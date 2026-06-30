import os
import sys

# Adiciona o diretório principal ao path para poder importar core e run_validation
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import run_validation as rv
import core.pipeline_common as common

def validate_single_round(
    round_id_target: int,
    settings_path: str = "config/settings.json",
    snippets_path: str = "candidates_subset.json",
    ledger_path: str = "checkpoint_ledger_phase2.csv",
    baseline_cache_path: str = "baseline_cache.json"
) -> None:
    # Garante que o SonarQube está online (se configurado)
    rv.ensure_sonarqube_online(settings_path)
    
    settings, _ = common.load_config(settings_path, snippets_path)
    checkout_dir = common.get_checkout_dir(settings)
    
    print(f"=== Iniciando Validacao para o Round {round_id_target} ===")
    print(f"Checkout Dir: {checkout_dir}")
    print(f"Ledger File: {ledger_path}")
    
    # Executa a validacao apenas para a rodada especificada
    rv.run_validation(
        settings_path=settings_path,
        snippets_path=snippets_path,
        ledger_path=ledger_path,
        baseline_cache_path=baseline_cache_path,
        checkout_dir_override=checkout_dir,
        only_round_ids={round_id_target},
        freeze_baseline=False
    )
    print("=== Validacao Concluida ===")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Uso: python scratch/validate_single.py <round_id> [ledger_path]")
        sys.exit(1)
        
    rid = int(sys.argv[1])
    lpath = sys.argv[2] if len(sys.argv) > 2 else "checkpoint_ledger_phase2.csv"
    validate_single_round(round_id_target=rid, ledger_path=lpath)
