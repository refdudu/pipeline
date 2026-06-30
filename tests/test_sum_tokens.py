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
