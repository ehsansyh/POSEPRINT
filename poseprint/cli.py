"""Command line entry point.

    poseprint complex.pdb
    poseprint poses/*.pdb --out results --format fingerprint
    poseprint complex.pdb --mol2 complex.mol2 --full

One row per ligand for the fingerprint form, or the full per-residue table when --full
is given. Writing to a directory keeps one CSV per complex, writing to a file collects
everything into it.
"""
from __future__ import annotations

import argparse
import glob
import sys
import warnings
from pathlib import Path

import pandas as pd

from . import __version__
from .core import INTERACTION_COLUMNS, analyse


def _expand(patterns):
    out = []
    for pattern in patterns:
        hits = sorted(glob.glob(pattern))
        if hits:
            out.extend(Path(h) for h in hits)
        elif Path(pattern).is_file():
            out.append(Path(pattern))
        else:
            print("  no such file, %s" % pattern, file=sys.stderr)
    return out


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="poseprint",
        description="Protein-ligand interaction fingerprints from docked poses.")
    p.add_argument("complexes", nargs="+",
                   help="one or more complexes in PDB format, wildcards allowed")
    p.add_argument("--mol2", default=None,
                   help="the matching MOL2, only valid with a single complex. "
                        "Built with Open Babel when left out")
    p.add_argument("--out", "-o", default=None,
                   help="a CSV file, or a directory to fill with one CSV per complex. "
                        "Prints to the screen when left out")
    p.add_argument("--full", action="store_true",
                   help="write the whole per-residue table instead of the summary")
    p.add_argument("--all-residues", action="store_true",
                   help="with --full, keep residues that make no contact at all")
    p.add_argument("--diagram", "-d", nargs="?", const=".", default=None,
                   metavar="DIR",
                   help="also draw the two-dimensional interaction diagram into DIR, "
                        "or the working directory when DIR is left out")
    p.add_argument("--diagram-format", default="pdf",
                   choices=["pdf", "png", "svg", "eps"],
                   help="the file type for --diagram, pdf by default so it stays sharp "
                        "at any size in a manuscript")
    p.add_argument("--log", nargs="?", const=".", default=None, metavar="DIR",
                   help="also write the report for each complex into DIR, as the text "
                        "POSEPRINT writes and as a CSV of one row per contact")
    p.add_argument("--quiet", "-q", action="store_true", help="no progress output")
    p.add_argument("--version", action="version", version="POSEPRINT %s" % __version__)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    warnings.filterwarnings("ignore")
    paths = _expand(args.complexes)
    if not paths:
        print("  nothing to do", file=sys.stderr)
        return 1
    if args.mol2 and len(paths) > 1:
        print("  --mol2 only makes sense with a single complex", file=sys.stderr)
        return 1

    out = Path(args.out) if args.out else None
    per_file = out is not None and (out.is_dir() or out.suffix == "")
    if per_file:
        out.mkdir(parents=True, exist_ok=True)

    collected, failures = [], 0
    for i, path in enumerate(paths, 1):
        if not args.quiet:
            print("  [%d/%d] %s" % (i, len(paths), path.name), file=sys.stderr)
        try:
            results = analyse(path, mol2=args.mol2, verbose=False)
        except Exception as exc:                       # keep going over a batch
            failures += 1
            print("    failed, %s: %s" % (type(exc).__name__, exc), file=sys.stderr)
            continue
        if not results:
            print("    no ligand found in %s" % path.name, file=sys.stderr)
            continue
        if args.diagram is not None:
            target = Path(args.diagram)
            target.mkdir(parents=True, exist_ok=True)
            try:
                from .diagram import draw_diagram
                import matplotlib.pyplot as plt
                figure_path = target / ("%s_interactions.%s"
                                        % (path.stem, args.diagram_format))
                fig = draw_diagram(path, out=figure_path, mol2=args.mol2,
                                   title=path.stem)
                plt.close(fig)
                if not args.quiet:
                    print("    diagram -> %s" % figure_path, file=sys.stderr)
            except Exception as exc:
                print("    the diagram could not be drawn, %s: %s"
                      % (type(exc).__name__, exc), file=sys.stderr)

        if args.log is not None:
            target = Path(args.log)
            target.mkdir(parents=True, exist_ok=True)
            for res in results:
                stem = "%s_%s" % (path.stem, res.ligand)
                (target / (stem + ".log")).write_text(res.log, encoding="utf-8")
                res.log_table().to_csv(target / (stem + "_interactions.csv"),
                                       index=False)
            if not args.quiet:
                print("    reports -> %s" % target, file=sys.stderr)

        for res in results:
            if args.full:
                frame = res.table if args.all_residues else res.contact_table()
                frame = frame.copy()
            else:
                counts = res.counts()
                frame = pd.DataFrame([{
                    "complex": path.name,
                    "ligand": res.ligand,
                    "hbond_count": counts.get("H-bonds", 0),
                    "hbond_residues": ";".join(res.residues("H-bonds")),
                    "hydrophobic_count": counts.get("Hydrophobic interaction", 0),
                    "hydrophobic_residues": ";".join(
                        res.residues("Hydrophobic interaction")),
                    **{c.lower().replace(" ", "_"): counts.get(c, 0)
                       for c in INTERACTION_COLUMNS
                       if c not in ("H-bonds", "Hydrophobic interaction")},
                }])
            if "complex" not in frame.columns:
                frame.insert(0, "complex", path.name)
            if "ligand" not in frame.columns:
                frame.insert(1, "ligand", res.ligand)
            if per_file:
                target = out / ("%s_%s.csv" % (path.stem, res.ligand))
                frame.to_csv(target, index=False)
            else:
                collected.append(frame)

    if collected:
        table = pd.concat(collected, ignore_index=True)
        if out is not None:
            out.parent.mkdir(parents=True, exist_ok=True)
            table.to_csv(out, index=False)
            if not args.quiet:
                print("  wrote %s, %d row(s)" % (out, len(table)), file=sys.stderr)
        else:
            print(table.to_string(index=False))
    elif per_file and not args.quiet:
        print("  wrote one CSV per complex into %s" % out, file=sys.stderr)

    if failures and not args.quiet:
        print("  %d complex(es) failed" % failures, file=sys.stderr)
    return 0 if failures == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
