import os
import sys
import logging
import json
from typing import Any, Dict, List, Tuple, Optional

# Standard personas. Only the identity/role preamble varies across levels; the
# refactoring base instruction is identical for all (so the persona is the only
# manipulated variable). P-2 is intentionally empty (neutral reference).
PERSONAS = {
    "P-1": "You are a master artisan baker with decades of experience crafting sourdough breads and traditional pastries. Baking is your entire world. You have never studied or practiced software engineering, programming, or computer science, and you have no knowledge of them whatsoever. Stay fully in character as this baker.",
    "P-2": "",
    "P-3": "You are a software developer.",
    "P-4": "You are a Senior Java Software Architect and a recognized Clean Code expert, with deep mastery of object-oriented design and software maintainability."
}

_logger = None

def get_logger() -> logging.Logger:
    global _logger
    if _logger is None:
        os.makedirs("logs", exist_ok=True)
        _logger = logging.getLogger("pipeline")
        _logger.setLevel(logging.DEBUG)
        _logger.propagate = False
        if not _logger.handlers:
            file_handler = logging.FileHandler("logs/pipeline.log", encoding="utf-8")
            file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            file_handler.setFormatter(file_formatter)
            _logger.addHandler(file_handler)
            
            console_handler = logging.StreamHandler(sys.stdout)
            console_handler.setLevel(logging.INFO)
            console_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            console_handler.setFormatter(console_formatter)
            _logger.addHandler(console_handler)
    return _logger

def load_config(settings_path: str, snippets_path: str) -> Tuple[dict, list]:
    if not os.path.exists(settings_path):
        raise FileNotFoundError(f"Settings file not found: {settings_path}")
    if not os.path.exists(snippets_path):
        raise FileNotFoundError(f"Snippets file not found: {snippets_path}")

    with open(settings_path, "r", encoding="utf-8") as f:
        settings = json.load(f)
    with open(snippets_path, "r", encoding="utf-8") as f:
        snippets = json.load(f)

    if isinstance(snippets, dict):
        flat_snippets = []
        for key in snippets:
            if isinstance(snippets[key], list):
                flat_snippets.extend(snippets[key])
        snippets = flat_snippets

    return settings, snippets

def resolve_personas(personas_override: Optional[Dict[str, str]]) -> dict:
    return personas_override if personas_override is not None else PERSONAS

def resolve_replicas(replicas_override: Optional[int]) -> List[int]:
    replicas_limit = 2 if replicas_override is None else replicas_override
    return list(range(1, replicas_limit + 1))

def build_rounds(snippets: list, personas: dict, replicas: List[int]) -> list:
    rounds_metadata = []
    id_rodada = 1
    for snippet in snippets:
        for persona_key, persona_preamble in personas.items():
            for replica in replicas:
                rounds_metadata.append({
                    "id_rodada": id_rodada,
                    "snippet": snippet,
                    "persona_key": persona_key,
                    "persona_preamble": persona_preamble,
                    "replica": replica
                })
                id_rodada += 1
    return rounds_metadata

def get_checkout_dir(settings: dict) -> str:
    return os.path.join(settings["workspace_root"], "test_checkout")
