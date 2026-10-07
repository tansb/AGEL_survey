"""
batch_runner.py — process the next N pending targets for one proposal, end to end.

For each target: stage raw data (reuse HST_data cache if present, else fall back to
a real MAST download), run hst_reduction steps 2-4 (drizzle/reproject/green), then
hst_products steps 1+3 (cutouts + postage stamps; offset-fitting/step 2 skipped by
design for this batch run). Progress is tracked in batch_processing_tracker.csv,
updated after every target so partial progress is always visible.

main() only ever picks up tracker rows already marked 'pending' -- it does not
create them. Run --seed first (once per proposal, and again any time targets
are added to an in-progress proposal's CSV, e.g. an ongoing GAP program) to add
any catalogue_objname entries missing from the tracker as new 'pending' rows.

Usage:
    python batch_runner.py --proposal 17307 --seed          # add new targets as 'pending'
    python batch_runner.py --proposal 17307 --batch-size 10 # process the next 10 pending
"""

import argparse
import shutil
import traceback
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

import hst_reduction
import hst_products
from config import HST_DATA_CACHE, MAIN_DIR, POSTAGE_STAMP_CONFIG_CSV

TRACKER_PATH = Path(__file__).parent / 'batch_processing_tracker.csv'

PROPOSAL_META = {
    '15867': {'filters': ['F140W'],          'camera': 'WFC3'},
    '16773': {'filters': ['F140W', 'F200LP'], 'camera': 'WFC3'},
    '17307': {'filters': ['F606W'],          'camera': 'ACS'},
}

PROPOSAL_CSVS = {
    '15867': MAIN_DIR / '15867_targets.csv',
    '16773': MAIN_DIR / '16773_targets.csv',
    '17307': MAIN_DIR / '17307_targets.csv',
}


def _now():
    return datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')


def _band_glob(band):
    return '*flt.fits' if band == 'F140W' else '*flc.fits'


def stage_raw_data(objname, catalogue_objname, proposal_id, band):
    """Copy cached raw files from HST_DATA_CACHE if present; else download from MAST.
    Returns 'cached' or 'downloaded'."""
    dest_dir = MAIN_DIR / f'{objname}/HST/{proposal_id}_{band}/raw_data'
    src_dir = HST_DATA_CACHE / objname / 'HST' / f'{proposal_id}_{band}' / 'raw_data'

    cached_files = sorted(src_dir.glob(_band_glob(band))) if src_dir.exists() else []
    if cached_files:
        dest_dir.mkdir(parents=True, exist_ok=True)
        for f in cached_files:
            shutil.copy(str(f), str(dest_dir))
        return 'cached'

    hst_reduction.download(objname, catalogue_objname, band, proposal_id)
    return 'downloaded'


def load_alt_run_targets():
    """objnames with alt_run=True in postage_stamp_config.csv — their postage stamps
    are built from *_L3.fits / *_green_L3.fits files produced by a separate,
    standalone notebook workflow rather than this pipeline's own drizzle output."""
    ps_cfg = pd.read_csv(POSTAGE_STAMP_CONFIG_CSV)
    return set(ps_cfg.loc[ps_cfg['alt_run'] == True, 'objname'])


def stage_alt_run_files(objname, proposal_id, filters):
    """Copy the pre-existing alt_run L3 files from HST_DATA_CACHE into HST_workingdata
    so hst_products.run_step3's alt_run branch finds what it expects. filters[0] is
    treated as the red/IR band (paired with the green file), matching how
    _make_color_postage's alt_run branch reads them."""
    copied = []
    for i, band in enumerate(filters):
        src_dir = HST_DATA_CACHE / objname / 'HST' / f'{proposal_id}_{band}'
        dest_dir = MAIN_DIR / f'{objname}/HST/{proposal_id}_{band}'
        dest_dir.mkdir(parents=True, exist_ok=True)

        names = [f'{objname}_{band}_L3.fits']
        if i == 0:
            names.append(f'{objname}_{band}_green_L3.fits')

        for name in names:
            src = src_dir / name
            if src.exists():
                shutil.copy(str(src), str(dest_dir / name))
                copied.append(name)
    return copied


def ensure_postage_stamp_config_row(objname):
    """Every target should have a postage_stamp_config.csv row so its thumb_size/clims
    are individually tunable later, even though hst_products.run_step3 already falls
    back to defaults (thumb_size=20, clims_set=1) for targets missing from the CSV.
    No-op if the target already has a row."""
    ps_cfg = pd.read_csv(POSTAGE_STAMP_CONFIG_CSV)
    if objname in set(ps_cfg['objname']):
        return False

    new_row = {col: '' for col in ps_cfg.columns}
    new_row.update({
        'objname': objname,
        'thumb_size_arcsec': 20,
        'clims_set': 1,
        'skip': False,
        'alt_run': False,
    })
    ps_cfg = pd.concat([ps_cfg, pd.DataFrame([new_row])], ignore_index=True)
    ps_cfg.to_csv(POSTAGE_STAMP_CONFIG_CSV, index=False)
    return True


def process_target(row, alt_run_targets):
    objname = row['objname']
    catalogue_objname = row['catalogue_objname']
    proposal_id = row['proposal_id']
    filters = row['filters'].split(',')
    camera = row['camera']

    notes = []

    # 1. Stage raw data per band (cache-first, MAST fallback)
    for band in filters:
        source = stage_raw_data(objname, catalogue_objname, proposal_id, band)
        notes.append(f'{band}:{source}')

    # 1b. alt_run targets also need pre-existing L3 files for the postage stamp step
    if objname in alt_run_targets:
        copied = stage_alt_run_files(objname, proposal_id, filters)
        notes.append(f'alt_run_files:{len(copied)}')

    # 2. Drizzle (reduction step 2) — only element [0] (objname) is actually used
    # internally by run_step2/ir_drizzler/uv_drizzler.
    target_tuple = (objname, catalogue_objname, None, None, None, None)
    hst_reduction.run_step2([target_tuple], filters, proposal_id, camera, interactive=False)

    # 3. Reproject + 4. Green image (scoped to this target only)
    hst_reduction.run_step3(target_filter=objname)
    hst_reduction.run_step4(target_filter=objname)

    # 5. Cutouts (products step 1) — needs the full original CSV row as a dict
    hst_products.run_step1([row['_csv_row']], proposal_id, filters, camera)

    # 6. Postage stamp (products step 3), offset-fitting (step 2) skipped by design
    hst_products.run_step3(target_filter=objname)

    # Register in postage_stamp_config.csv (default thumb_size/clims) if not already
    # present, so it shows up for later per-target tuning.
    if ensure_postage_stamp_config_row(objname):
        notes.append('added_to_ps_config')

    return ','.join(notes)


def seed_tracker(proposal_id):
    """Add any targets from PROPOSAL_CSVS[proposal_id] missing from the tracker
    as new 'pending' rows. Idempotent — safe to run again after the catalogue
    CSV gains targets (e.g. an ongoing GAP program like 17307), since it only
    adds rows for (proposal_id, objname) pairs not already present; it never
    touches an existing row's status. 'B' entries are skipped, same convention
    as hst_reduction.run_step1/run_step2 (they're never independently
    processed — see README's 'Required input files' note on B entries).

    Returns the number of rows added.
    """
    tracker = pd.read_csv(TRACKER_PATH, dtype=str).fillna('') if TRACKER_PATH.exists() \
        else pd.DataFrame(columns=['proposal_id', 'objname', 'catalogue_objname', 'filters',
                                    'camera', 'batch_number', 'status', 'raw_data_source',
                                    'started_at', 'completed_at', 'notes'])
    existing = set(zip(tracker['proposal_id'], tracker['objname'])) if len(tracker) else set()

    meta = PROPOSAL_META[proposal_id]
    csv_df = pd.read_csv(PROPOSAL_CSVS[proposal_id])

    new_rows = []
    for _, r in csv_df.iterrows():
        objname = r['objname']
        if objname.endswith('B'):
            continue
        if (proposal_id, objname) in existing:
            continue
        new_rows.append({
            'proposal_id': proposal_id, 'objname': objname,
            'catalogue_objname': r['catalogue_objname'],
            'filters': ','.join(meta['filters']), 'camera': meta['camera'],
            'batch_number': '', 'status': 'pending', 'raw_data_source': '',
            'started_at': '', 'completed_at': '', 'notes': '',
        })

    if new_rows:
        tracker = pd.concat([tracker, pd.DataFrame(new_rows)], ignore_index=True)
        tracker.to_csv(TRACKER_PATH, index=False)
    return len(new_rows)


def main():
    parser = argparse.ArgumentParser(description='Process next N pending targets for one proposal')
    parser.add_argument('--proposal', required=True, choices=list(PROPOSAL_META))
    parser.add_argument('--batch-size', type=int, default=10)
    parser.add_argument('--seed', action='store_true',
                         help="Add new targets from the proposal's CSV to the tracker as "
                              "'pending' rows, then exit without processing anything.")
    args = parser.parse_args()

    if args.seed:
        n = seed_tracker(args.proposal)
        print(f"Added {n} new pending row(s) for proposal {args.proposal}."
              if n else f"No new targets for proposal {args.proposal} — tracker already up to date.")
        return

    tracker = pd.read_csv(TRACKER_PATH, dtype=str).fillna('')

    existing_batches = tracker.loc[tracker['batch_number'] != '', 'batch_number']
    # astype(float) first: tolerates "8" and "8.0" alike, in case the tracker was ever
    # re-saved without dtype=str (pandas then infers batch_number as float64).
    batch_number = (int(existing_batches.astype(float).max()) + 1) if len(existing_batches) else 1

    pending_mask = (tracker['proposal_id'] == args.proposal) & (tracker['status'] == 'pending')
    pending_idx = tracker[pending_mask].index[:args.batch_size]

    if len(pending_idx) == 0:
        print(f"No pending targets left for proposal {args.proposal}.")
        return

    print(f"=== Batch {batch_number}: {len(pending_idx)} target(s) from proposal {args.proposal} ===")

    csv_df = pd.read_csv(PROPOSAL_CSVS[args.proposal])
    csv_by_name = {r['objname']: r for _, r in csv_df.iterrows()}
    alt_run_targets = load_alt_run_targets()

    for idx in pending_idx:
        row = tracker.loc[idx].to_dict()
        objname = row['objname']
        print(f"\n--- {objname} ({args.proposal}) ---")

        tracker.loc[idx, 'batch_number'] = str(batch_number)
        tracker.loc[idx, 'status'] = 'in_progress'
        tracker.loc[idx, 'started_at'] = _now()
        tracker.to_csv(TRACKER_PATH, index=False)

        row['_csv_row'] = csv_by_name[objname].to_dict()

        try:
            notes = process_target(row, alt_run_targets)
            tracker.loc[idx, 'status'] = 'complete'
            tracker.loc[idx, 'raw_data_source'] = notes
            tracker.loc[idx, 'notes'] = ''
            print(f"  {objname}: complete ({notes})")
        except Exception as exc:
            tracker.loc[idx, 'status'] = 'failed'
            tracker.loc[idx, 'notes'] = f'{type(exc).__name__}: {exc}'
            print(f"  {objname}: FAILED — {type(exc).__name__}: {exc}")
            traceback.print_exc()
        finally:
            tracker.loc[idx, 'completed_at'] = _now()
            tracker.to_csv(TRACKER_PATH, index=False)

    print(f"\n=== Batch {batch_number} complete ===")
    summary = tracker[tracker['batch_number'] == str(batch_number)]['status'].value_counts()
    print(summary.to_string())


if __name__ == '__main__':
    main()
