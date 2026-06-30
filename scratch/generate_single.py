import os
import sys
import json
import logging
from typing import Dict, Any, Optional

# Adiciona o diretório principal ao path para poder importar core
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import core.pipeline_common as common
from generate_refactorings import init_gemini_clients, RefactorResult
from google.genai import types

def generate_single_round(
    round_id_target: int,
    settings_path: str = "config/settings.json",
    snippets_path: str = "candidates_subset.json",
    force_java7: bool = True
) -> None:
    logger = common.get_logger()
    logger.setLevel(logging.INFO)
    
    settings, snippets = common.load_config(settings_path, snippets_path)
    personas = common.resolve_personas(None)
    replicas = common.resolve_replicas(None)
    rounds_metadata = common.build_rounds(snippets, personas, replicas)
    
    target_round = None
    for r in rounds_metadata:
        if r["id_rodada"] == round_id_target:
            target_round = r
            break
            
    if not target_round:
        logger.error(f"Erro: Round ID {round_id_target} nao encontrado em {snippets_path}.")
        sys.exit(1)
        
    snippet = target_round["snippet"]
    persona_key = target_round["persona_key"]
    persona_preamble = target_round["persona_preamble"]
    replica = target_round["replica"]
    snippet_id = snippet["trecho"]
    
    logger.info(f"=== Iniciando Geracao para o Round {round_id_target} ===")
    logger.info(f"Snippet: {snippet_id}")
    logger.info(f"Persona: {persona_key}")
    logger.info(f"Replica: {replica}")
    
    checkout_dir = common.get_checkout_dir(settings)
    refac_dir = os.path.join("refactored_code", snippet_id, persona_key, f"replica_{replica}")
    os.makedirs(refac_dir, exist_ok=True)
    
    refac_file = os.path.join(refac_dir, "code.java")
    response_file = os.path.join(refac_dir, "response.json")
    tokens_file = os.path.join(refac_dir, "tokens.json")
    
    rel_path = snippet["file_path"].replace("temp_workspace/finder_Math_1/", "")
    source_file_path = os.path.join(checkout_dir, rel_path)
    
    if not os.path.exists(source_file_path):
        logger.error(f"Arquivo original nao encontrado em {source_file_path}.")
        sys.exit(1)
        
    with open(source_file_path, "r", encoding="utf-8") as f:
        original_code = f.read()
        
    # Salva original.java se nao existir
    original_file = os.path.join("refactored_code", snippet_id, "original.java")
    if not os.path.exists(original_file):
        with open(original_file, "w", encoding="utf-8") as f_orig:
            f_orig.write(original_code)
            
    ref_type = snippet.get("refatoracao_tipo", "Refactor")
    refactor_instruction = f"Refactor the following Java class to improve its readability and maintainability without changing its observable behavior. Apply the '{ref_type}' refactoring throughout the entire class, wherever it is appropriate, not only in a single location. The host project targets Java 5 (source/target 1.5); the refactored code MUST compile under Java 5 - do NOT use lambdas, method references (::), streams, the diamond operator (<>), try-with-resources, multi-catch, var, or any language feature newer than Java 5. Return the complete refactored class."
    formatting_instruction = "\nIMPORTANT: The refactored Java code in your response MUST be properly formatted with standard indentation and newlines. Do not compress the code into a single line."
    
    full_instruction = refactor_instruction + formatting_instruction

    clients = init_gemini_clients(settings)
    if not clients:
        logger.error("Nenhum cliente Gemini inicializado.")
        sys.exit(1)
        
    client = clients[0] # Usa o primeiro cliente disponivel
    
    logger.info("Chamando Gemini API (com retentativas para 503)...")
    import time
    from google.genai.errors import APIError
    
    response = None
    for attempt in range(1, 6):
        try:
            response = client.models.generate_content(
                model="gemini-3.5-flash",
                contents=f"{full_instruction}\n\n{original_code}",
                config=types.GenerateContentConfig(
                    response_mime_type="application/json",
                    response_schema=RefactorResult,
                    system_instruction=persona_preamble or None,
                )
            )
            break
        except Exception as e:
            is_transient = False
            if hasattr(e, "code") and e.code in (429, 503):
                is_transient = True
            elif "503" in str(e) or "429" in str(e) or "demand" in str(e).lower():
                is_transient = True
                
            if is_transient and attempt < 5:
                wait_sec = attempt * 5
                logger.warning(f"Tentativa {attempt} falhou ({e}). Aguardando {wait_sec}s antes de tentar novamente...")
                time.sleep(wait_sec)
            else:
                raise e
                
    # Salva response cru
    raw_json_str = response.model_dump_json(indent=2)
    with open(response_file, "w", encoding="utf-8") as f_raw:
        f_raw.write(raw_json_str)
        
    # Salva codigo refatorado
    result = response.parsed
    if result is None or not getattr(result, "refactored_code", None):
        logger.error("Erro: Resposta estruturada do Gemini vazia ou invalida.")
        sys.exit(1)
        
    refactored_code = result.refactored_code
    with open(refac_file, "w", encoding="utf-8") as f_out:
        f_out.write(refactored_code)
        
    # Salva tokens
    usage = response.usage_metadata
    tokens_data = {
        "prompt_tokens": usage.prompt_token_count if usage else 0,
        "candidates_tokens": usage.candidates_token_count if usage else 0,
        "total_tokens": usage.total_token_count if usage else 0
    }
    with open(tokens_file, "w", encoding="utf-8") as f_tok:
        json.dump(tokens_data, f_tok, indent=2)
        
    logger.info(f"Sucesso! Codigo refatorado salvo em: {refac_file}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Uso: python scratch/generate_single.py <round_id> [--no-java7]")
        sys.exit(1)
        
    rid = int(sys.argv[1])
    no_java7 = "--no-java7" in sys.argv
    generate_single_round(round_id_target=rid, force_java7=not no_java7)
