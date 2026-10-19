"""Command-line entry point:  ``kcwilensfit <command>``.

    kcwilensfit gui [config.json]       open the graphical interface
    kcwilensfit run config.json         run every stage headlessly and write the outputs
    kcwilensfit init-config out.json    write a config file with all defaults (and descriptions in docs)
    kcwilensfit demo [folder] [--run]   create synthetic compound-lens cubes + config (and run them)
    kcwilensfit fields                  print every setting with its default and explanation
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import textwrap


def _cmd_run(args) -> int:
    from .pipeline import Pipeline
    from .schema import load_config, validate
    cfg = load_config(args.config)
    if args.output:
        cfg["output_dir"] = args.output
    if cfg.pop("_unknown_keys", None):
        pass
    probs = validate(cfg)
    if probs:
        print("Configuration problems:\n  - " + "\n  - ".join(probs), file=sys.stderr)
        return 2
    pipe = Pipeline(cfg)
    pipe.run_all(save=not args.no_save)
    print("\n=== Summary ===")
    for k, r in pipe.results.items():
        extra = f", sigma={r.sigma:.0f}±{r.sigma_err:.0f} km/s" if r.sigma == r.sigma else ""
        print(f"{k:12s} {r.role:9s} {r.status:12s} z={r.z:.5f} ± {r.z_err:.5f} ({r.method}){extra}")
    for m in pipe.system_messages:
        print("  *", m)
    if not args.no_save:
        print(f"\nOutputs written to {os.path.abspath(cfg['output_dir'])}")
    return 0


def _cmd_gui(args) -> int:
    from .gui.app import launch
    launch(args.config)
    return 0


def _cmd_init(args) -> int:
    from .schema import defaults, save_config
    cfg = defaults()
    save_config(cfg, args.path)
    print(f"wrote {args.path}")
    return 0


def _cmd_fields(args) -> int:
    from .schema import FIELDS, STAGES
    names = dict(STAGES)
    stage = None
    for f in FIELDS:
        if f.stage != stage:
            stage = f.stage
            print(f"\n=== {names.get(stage, stage)} ===")
        print(f"\n{f.key}  [{f.kind}]  default: {f.default!r}")
        print(textwrap.indent(textwrap.fill(f.help, 96), "    "))
    return 0


def _cmd_demo(args) -> int:
    from .synthetic import make_demo_cubes
    from .schema import defaults, save_config
    folder = os.path.abspath(args.folder)
    d = make_demo_cubes(folder)
    cfg = defaults()
    cfg.update(blue_cube=d["paths"]["blue"], red_cube=d["paths"]["red"], output_dir=os.path.join(folder, "output"),
               target_name="synthetic_compound_lens", apertures=[dict(a) for a in d["apertures"]])
    from .schema import normalise_aperture
    cfg["apertures"] = [normalise_aperture(a) for a in cfg["apertures"]]
    cpath = os.path.join(folder, "demo_config.json")
    save_config(cfg, cpath)
    print(f"Synthetic cubes and config written to {folder}")
    print("Truth:", json.dumps(d["truth"]))
    if args.run:
        ns = argparse.Namespace(config=cpath, output=None, no_save=False)
        return _cmd_run(ns)
    print(f"Run it with:  kcwilensfit run {cpath}    or open it in the GUI:  kcwilensfit gui {cpath}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="kcwilensfit", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("gui", help="open the graphical interface")
    g.add_argument("config", nargs="?", default=None)
    g.set_defaults(fn=_cmd_gui)
    r = sub.add_parser("run", help="run all stages headlessly")
    r.add_argument("config")
    r.add_argument("-o", "--output", default=None, help="override the output directory")
    r.add_argument("--no-save", action="store_true")
    r.set_defaults(fn=_cmd_run)
    i = sub.add_parser("init-config", help="write a default config file")
    i.add_argument("path", nargs="?", default="kcwilensfit_config.json")
    i.set_defaults(fn=_cmd_init)
    f = sub.add_parser("fields", help="list every setting with its default and description")
    f.set_defaults(fn=_cmd_fields)
    d = sub.add_parser("demo", help="make synthetic compound-lens cubes (and optionally run them)")
    d.add_argument("folder", nargs="?", default="kcwilensfit_demo")
    d.add_argument("--run", action="store_true")
    d.set_defaults(fn=_cmd_demo)
    args = ap.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
