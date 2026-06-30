import os
import sys
import time
import json
import logging
from typing import Any, Dict, List, Optional

from google import genai
from google.genai import types
from pydantic import BaseModel, Field

import core.pipeline_common as common

class RefactorResult(BaseModel):
    refactored_code: str = Field(description="The refactored Java code. MUST preserve newlines and standard indentation.")

def _parse_env_file(path: str) -> Dict[str, str]:
    """Parse a simple KEY=VALUE .env file into a dict (ignores blanks/comments)."""
    values: Dict[str, str] = {}
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, v = line.split("=", 1)
            values[k.strip()] = v.strip().strip("'\" ")
    return values


def load_api_keys(settings: dict) -> List[str]:
    """Collect Gemini API keys for round-robin rotation. Supports multiple keys.
    Priority (first source with any key wins): settings 'api_keys' list >
    settings 'api_key' > env GEMINI_API_KEYS (comma-separated) > env
    GEMINI_API_KEY > .env file (GEMINI_API_KEYS then GEMINI_API_KEY)."""
    keys: List[str] = []

    raw = settings.get("api_keys")
    if isinstance(raw, list):
        keys = list(raw)
    elif isinstance(raw, str) and raw.strip():
        keys = [raw]

    if not keys and settings.get("api_key"):
        keys = [settings["api_key"]]

    if not keys:
        if os.environ.get("GEMINI_API_KEYS"):
            keys = os.environ["GEMINI_API_KEYS"].split(",")
        elif os.environ.get("GEMINI_API_KEY"):
            keys = [os.environ["GEMINI_API_KEY"]]

    if not keys and os.path.exists(".env"):
        env_vals = _parse_env_file(".env")
        if env_vals.get("GEMINI_API_KEYS"):
            keys = env_vals["GEMINI_API_KEYS"].split(",")
        elif env_vals.get("GEMINI_API_KEY"):
            keys = [env_vals["GEMINI_API_KEY"]]

    # Normalize: strip whitespace/quotes, drop empties, dedupe preserving order.
    seen = set()
    clean: List[str] = []
    for k in keys:
        k = k.strip().strip("'\"")
        if k and k not in seen:
            seen.add(k)
            clean.append(k)
    return clean


def init_gemini_clients(settings: dict) -> List[genai.Client]:
    """Build one genai.Client per available API key (for round-robin rotation)."""
    keys = load_api_keys(settings)
    logger = common.get_logger()
    if not keys:
        logger.error("Error: no Gemini API key found (set GEMINI_API_KEYS or GEMINI_API_KEY).")
        sys.exit(1)
    logger.info(f"Loaded {len(keys)} Gemini API key(s) for round-robin rotation.")
    return [genai.Client(api_key=k) for k in keys]


def _is_retryable_error(e: Exception) -> bool:
    """True for transient errors worth retrying with another key: rate-limit/quota
    (HTTP 429) and transient server errors (HTTP 500/502/503/504, e.g. UNAVAILABLE
    'high demand')."""
    if getattr(e, "code", None) in (429, 500, 502, 503, 504):
        return True
    msg = str(e).lower()
    keywords = ("resource_exhausted", "rate limit", "rate-limit", "quota",
                "unavailable", "overloaded", "high demand", "try again",
                "503", "500", "502", "504")
    return any(k in msg for k in keywords)

def run_generation(
    settings_path: str = "config/settings.json",
    snippets_path: str = "candidates_subset.json",
    personas_override: Optional[Dict[str, str]] = None,
    replicas_override: Optional[int] = None
) -> None:
    logger = common.get_logger()
    logger.info("Starting Phase 1: Generation.")

    settings, snippets = common.load_config(settings_path, snippets_path)
    personas = common.resolve_personas(personas_override)
    replicas = common.resolve_replicas(replicas_override)
    rounds_metadata = common.build_rounds(snippets, personas, replicas)

    total_rounds = len(rounds_metadata)
    os.makedirs("refactored_code", exist_ok=True)

    clients = init_gemini_clients(settings)
    rotation_index = 0
    checkout_dir = common.get_checkout_dir(settings)

    for round_info in rounds_metadata:
        rid = round_info["id_rodada"]
        snippet = round_info["snippet"]
        persona_key = round_info["persona_key"]
        persona_preamble = round_info["persona_preamble"]
        replica = round_info["replica"]

        snippet_id = snippet["trecho"]
        refac_dir = os.path.join("refactored_code", snippet_id, persona_key, f"replica_{replica}")
        os.makedirs(refac_dir, exist_ok=True)

        refac_file = os.path.join(refac_dir, "code.java")
        response_file = os.path.join(refac_dir, "response.json")
        tokens_file = os.path.join(refac_dir, "tokens.json")

        if os.path.exists(refac_file) and os.path.exists(response_file):
            logger.info(f"Refactoring for Round {rid} already generated. Skipping API call.")
            continue

        logger.info(f"Generating Refactoring for Round {rid}/{total_rounds} (Snippet: {snippet_id}, Persona: {persona_key}, Replica: {replica})")

        try:
            rel_path = snippet["file_path"].replace("temp_workspace/finder_Math_1/", "")
            source_file_path = os.path.join(checkout_dir, rel_path)

            if not os.path.exists(source_file_path):
                raise FileNotFoundError(f"Source file not found at {source_file_path}. Run defects4j checkout in {checkout_dir} first!")

            with open(source_file_path, "r", encoding="utf-8") as f:
                original_code = f.read()

            # Save the original source once at the snippet root (shared by all
            # personas/replicas), e.g. refactored_code/<snippet>/original.java
            original_file = os.path.join("refactored_code", snippet_id, "original.java")
            if not os.path.exists(original_file):
                with open(original_file, "w", encoding="utf-8") as f_orig:
                    f_orig.write(original_code)

            ref_type = snippet.get("refatoracao_tipo", "Refactor")
            refactor_instruction = f"Refactor the following Java class to improve its readability and maintainability without changing its observable behavior. Apply the '{ref_type}' refactoring throughout the entire class, wherever it is appropriate, not only in a single location. Return the complete refactored class."
            formatting_instruction = "\nIMPORTANT: The refactored Java code in your response MUST be properly formatted with standard indentation and newlines. Do not compress the code into a single line."
            full_instruction = refactor_instruction + formatting_instruction
            print(persona_preamble)
            # Round-robin across API keys, with failover to the next key on a
            # rate-limit error. Only the round fails (.error) if every key is limited.
            n_keys = len(clients)
            base = rotation_index
            rotation_index = (rotation_index + 1) % n_keys
            response = None
            for attempt in range(n_keys):
                key_idx = (base + attempt) % n_keys
                try:
                    response = clients[key_idx].models.generate_content(
                        model="gemini-3.5-flash",
                        contents=f"{full_instruction}\n\n{original_code}",
                        config=types.GenerateContentConfig(
                            response_mime_type="application/json",
                            response_schema=RefactorResult,
                            system_instruction=persona_preamble or None,
                        )
                    )
                    break
                except Exception as api_exc:
                    if _is_retryable_error(api_exc) and attempt < n_keys - 1:
                        logger.warning(f"Round {rid}: key #{key_idx + 1}/{n_keys} failed ({api_exc}); trying next key.")
                        continue
                    raise

            # Save raw response
            raw_json_str = response.model_dump_json(indent=2)
            with open(response_file, "w", encoding="utf-8") as f_raw:
                f_raw.write(raw_json_str)

            # Parse and save code
            result = response.parsed
            if result is None or not getattr(result, "refactored_code", None):
                raise ValueError("Parsed structured response is empty or invalid.")

            refactored_code = result.refactored_code
            with open(refac_file, "w", encoding="utf-8") as f_out:
                f_out.write(refactored_code)

            # Save tokens
            usage = response.usage_metadata
            tokens_data = {
                "prompt_tokens": usage.prompt_token_count if usage else 0,
                "candidates_tokens": usage.candidates_token_count if usage else 0,
                "total_tokens": usage.total_token_count if usage else 0
            }
            with open(tokens_file, "w", encoding="utf-8") as f_tok:
                json.dump(tokens_data, f_tok, indent=2)

            logger.info(f"Saved refactored code, response and tokens to: {refac_dir}")

            # Clean up error file if it exists from a previous run
            err_file = refac_file + ".error"
            if os.path.exists(err_file):
                try:
                    os.remove(err_file)
                except Exception:
                    pass

        except Exception as e:
            logger.error(f"Error generating refactoring for Round {rid}: {e}")
            with open(refac_file + ".error", "w", encoding="utf-8") as f_err:
                f_err.write(str(e))

        print("sleep(60)")
        time.sleep(60)

if __name__ == "__main__":
    settings_arg = sys.argv[1] if len(sys.argv) > 1 else "config/settings.json"
    snippets_arg = sys.argv[2] if len(sys.argv) > 2 else "candidates_subset.json"
    run_generation(settings_path=settings_arg, snippets_path=snippets_arg)
