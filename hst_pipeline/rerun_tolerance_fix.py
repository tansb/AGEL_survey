"""
One-off re-run of the 18 batch_processing_tracker.csv targets that failed with
the Defect-4 "TweakReg produced NaN shift/rotation/scale" error, after adding
UV_DRIZZLE_DEFAULTS['tolerance'] = 3.0 to config.py (drizzlepac's default
tolerance was rejecting the final source-match step even when the coarse
2dhist shift search found hundreds of good matches). Reuses batch_runner's
own process_target() so tracking/notes stay consistent with the normal batch
flow, but targets this specific list instead of the next-N-pending queue.

Usage:
    python rerun_tolerance_fix.py
"""

import pandas as pd

import batch_runner
from batch_runner import TRACKER_PATH, PROPOSAL_CSVS, process_target, _now

TARGETS = [
    ('16773', 'AGEL010158-491738A'),
    ('16773', 'AGEL010238+015857A'),
    ('16773', 'AGEL024303-000600A'),
    ('16773', 'AGEL030022-500129A'),
    ('17307', 'AGEL011759-052718A'),
    ('17307', 'AGEL023211+001339A'),
    ('17307', 'AGEL053349-253654A'),
    ('17307', 'AGEL101103-001649A'),
    ('17307', 'AGEL111800-153227A'),
    ('17307', 'AGEL112354+505149A'),
    ('17307', 'AGEL144149+144121A'),
    ('17307', 'AGEL224504-501725A'),
    ('17307', 'AGEL233607-535236A'),
    ('17307', 'AGEL014433-114212A'),
    ('17307', 'AGEL012453-144303A'),
    ('17307', 'AGEL014919-134905A'),
    ('17307', 'AGEL064100+642915A'),
    ('17307', 'AGEL033717-315214A'),
]


def main():
    tracker = pd.read_csv(TRACKER_PATH, dtype=str).fillna('')
    existing_batches = tracker.loc[tracker['batch_number'] != '', 'batch_number']
    batch_number = (int(existing_batches.astype(float).max()) + 1) if len(existing_batches) else 1

    csv_cache = {}
    alt_run_targets = batch_runner.load_alt_run_targets()

    results = []
    for proposal_id, objname in TARGETS:
        matches = tracker.index[(tracker['proposal_id'] == proposal_id) & (tracker['objname'] == objname)]
        if len(matches) != 1:
            print(f"!! {objname} ({proposal_id}): expected exactly 1 tracker row, found {len(matches)} -- skipping")
            continue
        idx = matches[0]

        if proposal_id not in csv_cache:
            csv_df = pd.read_csv(PROPOSAL_CSVS[proposal_id])
            csv_cache[proposal_id] = {r['objname']: r for _, r in csv_df.iterrows()}

        row = tracker.loc[idx].to_dict()
        row['_csv_row'] = csv_cache[proposal_id][objname].to_dict()

        print(f"\n--- {objname} ({proposal_id}) ---")
        tracker.loc[idx, 'batch_number'] = str(batch_number)
        tracker.loc[idx, 'status'] = 'in_progress'
        tracker.loc[idx, 'started_at'] = _now()
        tracker.to_csv(TRACKER_PATH, index=False)

        try:
            notes = process_target(row, alt_run_targets)
            tracker.loc[idx, 'status'] = 'complete'
            tracker.loc[idx, 'raw_data_source'] = notes
            tracker.loc[idx, 'notes'] = ''
            print(f"  {objname}: complete ({notes})")
            results.append((objname, proposal_id, 'complete', ''))
        except Exception as exc:
            import traceback
            tracker.loc[idx, 'status'] = 'failed'
            tracker.loc[idx, 'notes'] = f'{type(exc).__name__}: {exc}'
            print(f"  {objname}: FAILED — {type(exc).__name__}: {exc}")
            traceback.print_exc()
            results.append((objname, proposal_id, 'failed', f'{type(exc).__name__}: {exc}'))
        finally:
            tracker.loc[idx, 'completed_at'] = _now()
            tracker.to_csv(TRACKER_PATH, index=False)

    print(f"\n=== Re-run batch {batch_number} complete ===")
    for objname, proposal_id, status, note in results:
        print(f"  {status:10s} {proposal_id} {objname} {note}")


if __name__ == '__main__':
    main()
