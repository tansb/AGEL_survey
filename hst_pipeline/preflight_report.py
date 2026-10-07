"""
preflight_report.py — Defect 1 (WCSNAME/HDRLET) preflight check across every
target in a proposal, using pipeline_qc.select_astrometric_reference().

For each target with raw exposures on disk for BOTH filters, reports which
filter carries the genuine a posteriori astrometric fit. Does not guess for
targets with no local data, or targets where the two filters' solution
classes are ambiguous (both fit / neither fit) -- those are flagged
separately so they aren't mistaken for "checked and fine".

Usage:
    python preflight_report.py --main-dir /path/to/MAIN_DIR \
        --targets-csv /path/to/16773_targets.csv \
        --proposal 16773 --ir-filter F140W --uv-filter F200LP
"""
import argparse
import csv
from pathlib import Path

import pipeline_qc


def check_target(main_dir, target, proposal, ir_filter, uv_filter):
    ir_dir = main_dir / target / 'HST' / f'{proposal}_{ir_filter}' / 'raw_data'
    uv_dir = main_dir / target / 'HST' / f'{proposal}_{uv_filter}' / 'raw_data'

    ir_raw = sorted(ir_dir.glob('*flt.fits')) if ir_dir.exists() else []
    uv_raw = sorted(uv_dir.glob('*flc.fits')) if uv_dir.exists() else []

    if not ir_raw and not uv_raw:
        return {'target': target, 'status': 'NO DATA', 'detail': f'{ir_dir} and {uv_dir} not found'}
    if not ir_raw:
        return {'target': target, 'status': 'NO DATA (IR only missing)', 'detail': f'{ir_dir} not found'}
    if not uv_raw:
        return {'target': target, 'status': 'NO DATA (UV only missing)', 'detail': f'{uv_dir} not found'}

    ir_check = pipeline_qc.classify_wcs(ir_raw[0])
    uv_check = pipeline_qc.classify_wcs(uv_raw[0])

    try:
        reference, _report = pipeline_qc.select_astrometric_reference(ir_raw, uv_raw)
        status = f'OK -- reference = {reference}'
    except pipeline_qc.AmbiguousAstrometricReference as exc:
        reference = 'AMBIGUOUS'
        status = 'AMBIGUOUS -- needs manual review'

    return {
        'target': target,
        'status': status,
        'ir_wcstype': ir_check['wcstype'],
        'ir_fit': ir_check['is_a_posteriori'],
        'uv_wcstype': uv_check['wcstype'],
        'uv_fit': uv_check['is_a_posteriori'],
        'reference': reference,
        'differs': ir_check['is_a_posteriori'] != uv_check['is_a_posteriori'],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--main-dir', required=True)
    ap.add_argument('--targets-csv', required=True)
    ap.add_argument('--proposal', default='16773')
    ap.add_argument('--ir-filter', default='F140W')
    ap.add_argument('--uv-filter', default='F200LP')
    args = ap.parse_args()

    main_dir = Path(args.main_dir)
    with open(args.targets_csv, encoding='utf-8-sig') as f:
        targets = [row['objname'] for row in csv.DictReader(f) if row.get('objname')]

    results = []
    for target in targets:
        if target.endswith('B'):
            continue
        results.append(check_target(main_dir, target, args.proposal, args.ir_filter, args.uv_filter))

    checked = [r for r in results if 'reference' in r]
    no_data = [r for r in results if 'reference' not in r]
    differs = [r for r in checked if r.get('differs')]

    print(f"{'target':<24} {'status':<28} {'IR (' + args.ir_filter + ')':<45} {'UV (' + args.uv_filter + ')':<45} differs?")
    print('-' * 160)
    for r in results:
        if 'reference' in r:
            print(f"{r['target']:<24} {r['status']:<28} "
                  f"{('FIT: ' + r['ir_wcstype']) if r['ir_fit'] else ('no fit: ' + r['ir_wcstype']):<45} "
                  f"{('FIT: ' + r['uv_wcstype']) if r['uv_fit'] else ('no fit: ' + r['uv_wcstype']):<45} "
                  f"{'*** YES ***' if r['differs'] else 'no'}")
        else:
            print(f"{r['target']:<24} {r['status']:<28} {r['detail']}")

    print()
    print(f"=== SUMMARY: {len(results)} targets in CSV (excluding B entries) ===")
    print(f"  Checked (both filters' raw data present locally): {len(checked)}")
    print(f"  No local data (need to check /Volumes/AGEL or wherever MAIN_DIR "
          f"actually points): {len(no_data)}")
    print(f"  Solution classes DIFFER between filters (same misalignment risk "
          f"as DESJ0206-0114): {len(differs)}")
    if differs:
        print(f"    -> {', '.join(r['target'] for r in differs)}")
    ambiguous = [r for r in checked if r['reference'] == 'AMBIGUOUS']
    if ambiguous:
        print(f"  AMBIGUOUS (both or neither filter fit -- needs manual review): "
              f"{len(ambiguous)} -> {', '.join(r['target'] for r in ambiguous)}")


if __name__ == '__main__':
    main()
