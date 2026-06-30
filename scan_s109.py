"""Measure ReplaceMagicNumber refactoring success via the SonarQube rule
java:S109 ("Magic numbers should not be used").

The total `code_smells` metric is blind to magic-number removal here (S109 is not
in the default profile, and the count is noisy). This re-scans only the
ReplaceMagicNumber rounds (source-only rule → NO compile/test needed; behaviour
was already validated in Phase 2) with S109 active and records the delta of S109
violations per round.

Robustness:
- Each scan gets its OWN working directory (avoids the scanner ProjectLock
  ".scannerwork/.sonar_lock" conflict between serial scans).
- The CSV is written incrementally and the run is resumable (already-recorded
  rounds are skipped), so a mid-run failure never loses progress.

Prereq: a Java quality profile with java:S109 active must be the default
(profile "Sonar way S109" was created for this). Output: s109_magicnumber.csv.

Usage:  .venv/bin/python scan_s109.py [checkout_dir] [out_csv]
"""

import os
import sys
import csv
import json
import subprocess

import core.pipeline_common as common
from core.sonarqube_anal import SonarQubeManager

FINDER_PREFIX = "temp_workspace/finder_Math_1/"
PERSONAS = ["P-1", "P-2", "P-3", "P-4"]
REPS = [1, 2]
RULE = "java:S109"
FIELDS = ["trecho", "persona", "replica", "s109_base", "s109_pos", "s109_delta"]


def main(checkout_dir="temp_workspace/test_checkout_w1", out_csv="s109_magicnumber.csv"):
    settings, _ = common.load_config("config/settings.json", "candidates_subset.json")
    rmn = json.load(open("candidates_subset.json"))["ReplaceMagicNumber"]

    sm = SonarQubeManager(settings["sonar_url"], settings["sonar_token"],
                          os.path.abspath(settings["sonar_scanner_path"]))
    classes = os.path.abspath(os.path.join(checkout_dir, "target", "classes"))

    iso = os.path.abspath(os.path.join(settings["workspace_root"], ".scanner_s109"))
    os.environ["SONAR_USER_HOME"] = iso + "_home"
    os.makedirs(os.environ["SONAR_USER_HOME"], exist_ok=True)

    # Resume: collect rounds already recorded.
    done = set()
    if os.path.exists(out_csv):
        with open(out_csv, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                done.add((row["trecho"], row["persona"], row["replica"]))
    f_out = open(out_csv, "a", newline="", encoding="utf-8")
    writer = csv.DictWriter(f_out, fieldnames=FIELDS)
    if not done:
        writer.writeheader(); f_out.flush()

    def restore(rel):
        subprocess.run(["git", "-C", checkout_dir, "checkout", "--", rel], capture_output=True)

    def s109(key, src):
        # Unique working dir per scan → no ProjectLock collision between scans.
        sm.scan(key, key, src, classes, work_dir=f"{iso}_w_{key}")
        if not sm.wait_for_task(key):
            print(f"  WARN: CE task not SUCCESS for {key}", flush=True)
        n = len(sm.search_issues(key, [RULE]))
        sm.delete_project(key)
        return n

    for snip in rmn:
        sid = snip["trecho"]
        rel = snip["file_path"].replace(FINDER_PREFIX, "")
        src = os.path.abspath(os.path.join(checkout_dir, rel))
        todo = [(p, r) for p in PERSONAS for r in REPS if (sid, p, str(r)) not in done]
        if not todo:
            print(f"{sid}: já completo, pulando", flush=True)
            continue

        restore(rel)
        base = s109(f"s109base_{sid}".replace("-", ""), src)
        print(f"{sid}: baseline S109={base}", flush=True)

        for p, r in todo:
            code = f"refactored_code/{sid}/{p}/replica_{r}/code.java"
            if not os.path.exists(code):
                row = {"trecho": sid, "persona": p, "replica": r,
                       "s109_base": base, "s109_pos": "", "s109_delta": ""}
            else:
                with open(src, "w", encoding="utf-8") as fh:
                    fh.write(open(code, encoding="utf-8").read())
                pos = s109(f"s109_{sid}_{p}_{r}".replace("-", ""), src)
                row = {"trecho": sid, "persona": p, "replica": r,
                       "s109_base": base, "s109_pos": pos, "s109_delta": pos - base}
                restore(rel)
                print(f"  {p} r{r}: pos={pos} delta={pos - base}", flush=True)
            writer.writerow(row); f_out.flush()
        restore(rel)

    f_out.close()
    print("Done.")


if __name__ == "__main__":
    ck = sys.argv[1] if len(sys.argv) > 1 else "temp_workspace/test_checkout_w1"
    out = sys.argv[2] if len(sys.argv) > 2 else "s109_magicnumber.csv"
    main(ck, out)
