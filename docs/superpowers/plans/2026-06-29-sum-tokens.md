# Sum Tokens Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create a Python script at the root of the workspace that recursively finds all `response.json` files within `refactored_code`, parses their token counts from `usage_metadata`, and sums them.

**Architecture:** Use `pathlib` to find all `response.json` files, load each file as JSON, extract `prompt_token_count` (input) and `candidates_token_count` (output) from `usage_metadata`, and display the sum.

**Tech Stack:** Python 3 (standard library only).

---

### Task 1: Create sum_tokens.py

**Files:**
- Create: `sum_tokens.py`
- Test: `tests/test_sum_tokens.py`

- [ ] **Step 1: Write the failing test**

Write a unit test that mocks `response.json` files and checks if the summing logic works correctly.

```python
import json
import tempfile
from pathlib import Path
import pytest
from sum_tokens import parse_and_sum_tokens

def test_parse_and_sum_tokens():
    with tempfile.TemporaryDirectory() as tmpdir:
        base_path = Path(tmpdir)
        
        # Create dummy structure
        d1 = base_path / "dir1"
        d1.mkdir()
        f1 = d1 / "response.json"
        
        # Valid JSON with usage metadata
        data1 = {
            "usage_metadata": {
                "prompt_token_count": 100,
                "candidates_token_count": 50
            }
        }
        f1.write_text(json.dumps(data1))
        
        # Valid JSON without usage metadata (should be handled gracefully)
        d2 = base_path / "dir2"
        d2.mkdir()
        f2 = d2 / "response.json"
        data2 = {
            "some_other_data": {}
        }
        f2.write_text(json.dumps(data2))
        
        # Invalid JSON (should be handled gracefully)
        d3 = base_path / "dir3"
        d3.mkdir()
        f3 = d3 / "response.json"
        f3.write_text("invalid json")
        
        inputs, outputs = parse_and_sum_tokens(base_path)
        assert inputs == 100
        assert outputs == 50
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_sum_tokens.py`
Expected: FAIL with "ModuleNotFoundError: No module named 'sum_tokens'" or similar.

- [ ] **Step 3: Write minimal implementation**

Create `sum_tokens.py` with the implementation of `parse_and_sum_tokens` and a CLI execution block.

```python
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_sum_tokens.py`
Expected: PASS

- [ ] **Step 5: Run the script on actual data**

Run: `python sum_tokens.py`
Expected: Outputs actual token sums.

- [ ] **Step 6: Commit**

```bash
git add sum_tokens.py tests/test_sum_tokens.py
git commit -m "feat: add sum_tokens script and tests"
```
