"""Phase 2 parallel orchestrator.

Runs run_validation across N workers, each with its own copy of the Defects4J
checkout (so concurrent compile/test/source-overwrite never collide) and its own
ledger shard. Baselines are precomputed serially first (shared, frozen cache) so
workers never race to scan/write the same snippet baseline. Shards are merged
into the final ledger at the end.

Usage:
  .venv/bin/python run_validation_parallel.py [settings] [snippets] [final_ledger] [num_workers] [max_per_worker]
"""

import os
import csv
import sys
import shutil
import subprocess
import multiprocessing as mp
from typing import Optional, List

import core.pipeline_common as common
import core.ledger
import run_validation as rv

WORKER_CHECKOUT_FMT = "test_checkout_w{n}"


def _ensure_worker_checkout(workspace_root: str, base_checkout: str, n: int) -> str:
    """Create (or reuse) worker N's checkout copy. Copies preserve .git + target/
    so each worker can git-restore independently and skip recompiling deps."""
    dst = os.path.join(workspace_root, WORKER_CHECKOUT_FMT.format(n=n))
    if os.path.isdir(dst):
        # Reuse existing copy; discard any leftover edits.
        subprocess.run(["git", "-C", dst, "checkout", "--", "."], capture_output=True, stdin=subprocess.DEVNULL)
        return dst
    print(f"[orchestrator] copying {base_checkout} -> {dst} ...", flush=True)
    shutil.copytree(base_checkout, dst, symlinks=True)
    return dst


def _worker(settings_path, snippets_path, ledger_path, baseline_cache_path,
            checkout_dir, round_ids, tag):
    """Subprocess entrypoint: validate only this worker's slice of round ids."""
    rv.run_validation(
        settings_path=settings_path,
        snippets_path=snippets_path,
        ledger_path=ledger_path,
        baseline_cache_path=baseline_cache_path,
        checkout_dir_override=checkout_dir,
        only_round_ids=set(round_ids),
        freeze_baseline=True,
        worker_tag=tag,
    )


def _merge_shards(shard_paths: List[str], final_path: str) -> None:
    """Merge worker ledger shards into one ledger, sorted by id_rodada, last write
    per id winning (supports resume)."""
    rows = {}
    for sp in shard_paths:
        if not os.path.exists(sp):
            continue
        with open(sp, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                try:
                    rid = int(row["id_rodada"])
                except (KeyError, ValueError, TypeError):
                    continue
                rows[rid] = row
    with open(final_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=core.ledger.LedgerManager.FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        for rid in sorted(rows):
            writer.writerow(rows[rid])
    print(f"[orchestrator] merged {len(rows)} rows -> {final_path}", flush=True)
    # Status histogram.
    hist = {}
    for r in rows.values():
        hist[r.get("status", "?")] = hist.get(r.get("status", "?"), 0) + 1
    print(f"[orchestrator] status counts: {hist}", flush=True)


def run_parallel(
    settings_path: str = "config/settings.json",
    snippets_path: str = "candidates_subset.json",
    final_ledger_path: str = "checkpoint_ledger_phase2.csv",
    num_workers: int = 3,
    baseline_cache_path: str = "baseline_cache.json",
    max_per_worker: Optional[int] = None,
) -> None:
    settings, snippets = common.load_config(settings_path, snippets_path)
    workspace_root = settings["workspace_root"]
    base_checkout = common.get_checkout_dir(settings)
    if not os.path.isdir(base_checkout):
        print(f"[orchestrator] base checkout missing: {base_checkout}", flush=True)
        sys.exit(1)

    # 1. Precompute all snippet baselines serially (shared, frozen cache).
    print("[orchestrator] precomputing baselines (serial)...", flush=True)
    rv.precompute_baselines(settings_path, snippets_path, baseline_cache_path, checkout_dir_override=base_checkout)

    # 2. Build the full round list and assign ids round-robin across workers
    #    (balances persona/refactoring-type mix per worker).
    rounds = common.build_rounds(snippets, common.resolve_personas(None), common.resolve_replicas(None))
    assignments = [[] for _ in range(num_workers)]
    for idx, r in enumerate(rounds):
        assignments[idx % num_workers].append(r["id_rodada"])
    if max_per_worker is not None:
        assignments = [a[:max_per_worker] for a in assignments]

    # 3. Spawn one worker process per checkout copy.
    shard_paths = []
    procs = []
    ctx = mp.get_context("spawn")
    for n in range(num_workers):
        checkout = _ensure_worker_checkout(workspace_root, base_checkout, n + 1)
        shard = f"{os.path.splitext(final_ledger_path)[0]}_w{n + 1}.csv"
        shard_paths.append(shard)
        tag = f"w{n + 1}"
        print(f"[orchestrator] {tag}: {len(assignments[n])} rounds, checkout={checkout}, shard={shard}", flush=True)
        p = ctx.Process(
            target=_worker,
            args=(settings_path, snippets_path, shard, baseline_cache_path, checkout, assignments[n], tag),
        )
        p.start()
        procs.append(p)

    for p in procs:
        p.join()

    # 4. Merge shards into the final ledger.
    _merge_shards(shard_paths, final_ledger_path)
    print("[orchestrator] done.", flush=True)


if __name__ == "__main__":
    settings_arg = sys.argv[1] if len(sys.argv) > 1 else "config/settings.json"
    snippets_arg = sys.argv[2] if len(sys.argv) > 2 else "candidates_subset.json"
    ledger_arg = sys.argv[3] if len(sys.argv) > 3 else "checkpoint_ledger_phase2.csv"
    workers_arg = int(sys.argv[4]) if len(sys.argv) > 4 else 3
    max_pw_arg = int(sys.argv[5]) if len(sys.argv) > 5 else None
    run_parallel(
        settings_path=settings_arg,
        snippets_path=snippets_arg,
        final_ledger_path=ledger_arg,
        num_workers=workers_arg,
        max_per_worker=max_pw_arg,
    )
