import json
import sys
from pathlib import Path

def parse_and_sum_tokens(root_dir: Path):
    total_inputs = 0
    total_outputs = 0
    
    # Recursively find response.json files
    for filepath in root_dir.glob("**/response.json"):
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                content = json.load(f)
                
            metadata = content.get("usage_metadata")
            if metadata:
                total_inputs += metadata.get("prompt_token_count") or 0
                total_outputs += metadata.get("candidates_token_count") or 0
        except Exception as e:
            # Silently ignore read/parse errors for robustness, or print a warning
            pass
            
    return total_inputs, total_outputs

if __name__ == "__main__":
    refactored_dir = Path(__file__).parent / "refactored_code"
    if not refactored_dir.exists():
        print(f"Error: {refactored_dir} does not exist.")
        sys.exit(1)
        
    inputs, outputs = parse_and_sum_tokens(refactored_dir)
    print(f"Total Input (Prompt) Tokens:  {inputs}")
    print(f"Total Output (Candidate) Tokens: {outputs}")
    print(f"Total Tokens:                   {inputs + outputs}")
