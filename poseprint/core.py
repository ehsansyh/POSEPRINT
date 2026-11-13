"""The interface most users want.

Point `analyse` at a docked complex and it returns one table per ligand, with a row for
every residue and a column for every interaction type. `fingerprint` reduces that table
to the two counts and the two residue lists that a scoring or learning model usually
wants, which is the form the CLINCH labels were built in.

Each call works in its own scratch directory. The detectors write intermediate ligand
files while they run, and the original code sent every call to one fixed directory, so
two analyses running at the same time could overwrite each other's intermediates.
"""
from __future__ import annotations

import os
import re
import shutil
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from . import analysis
from .convert import pdb_to_mol2

__all__ = ["InteractionResult", "analyse", "fingerprint", "INTERACTION_COLUMNS"]

# the per-residue columns run_analyze fills, in the order the table carries them
INTERACTION_COLUMNS = [
    "H-bonds",
    "Salt_bridge",
    "Water_bridge",
    "Pi-Stacking (T-shape)",
    "Pi-Stacking(F-F)",
    "Pi-Stacking(parallel)",
    "Pi-cation",
    "Hydrophobic interaction",
    "Halogen bonds",
]


# every line of the report has the same shape, a label, the partners, then the
# measurements, and sometimes the atoms that made the contact
LOG_LINE = re.compile(
    r"^(?P<label>.+?) between (?P<partners>.+?): (?P<fields>.+?)"
    r"(?: --> Interaction between atoms (?P<atoms>\(.*\)))?\s*$")


@dataclass
class InteractionResult:
    """Everything POSEPRINT found for one ligand in one complex."""

    pose: str
    ligand: str
    table: pd.DataFrame = field(repr=False)
    log: str = field(default="", repr=False)

    def _ligand_mask(self):
        label = "%s %s" % (self.ligand[:3], self.ligand[3:])
        return self.table["Residues"].astype(str).str.strip() == label

    def counts(self) -> dict:
        """How many contacts of each type the ligand makes.

        The table carries one row per residue plus one row for the ligand itself, and
        every contact appears twice, once against the residue and once against the
        ligand. The ligand row is therefore the total, and summing the whole column
        would double count. A homodimer makes that worse, because the residue label
        carries no chain, so the two copies of a residue share one label and one count.
        """
        mask = self._ligand_mask()
        rows = self.table[mask] if mask.any() else self.table
        out = {}
        for col in INTERACTION_COLUMNS:
            if col in self.table.columns:
                out[col] = int(pd.to_numeric(rows[col], errors="coerce").fillna(0).sum())
        return out

    def residues(self, column: str) -> list:
        """The protein residues that carry at least one contact of the given type.

        The ligand row is left out, and a residue is listed once even when a homodimer
        presents two copies of it under the same label.
        """
        if column not in self.table.columns:
            return []
        hits = pd.to_numeric(self.table[column], errors="coerce").fillna(0) > 0
        names = self.table.loc[hits & ~self._ligand_mask(), "Residues"]
        seen, out = set(), []
        for r in names:
            r = str(r).strip()
            if r not in seen:
                seen.add(r)
                out.append(r)
        return out

    def contact_table(self) -> pd.DataFrame:
        """The rows that carry at least one contact of any type."""
        cols = [c for c in INTERACTION_COLUMNS if c in self.table.columns]
        if not cols:
            return self.table.iloc[0:0]
        total = sum(pd.to_numeric(self.table[c], errors="coerce").fillna(0) for c in cols)
        keep = self.table.loc[total > 0].copy()
        keep["Residues"] = keep["Residues"].astype(str).str.strip()
        return keep.drop_duplicates(subset="Residues").reset_index(drop=True)

    def log_table(self) -> pd.DataFrame:
        """The report, one row per contact instead of one block of text.

        The per-residue table says how many contacts a residue makes. This says what
        each one is, which atoms made it and how far apart they are, which is what the
        report file carries and what a reader usually wants to quote.
        """
        rows = []
        for raw in self.log.splitlines():
            line = raw.strip()
            if not line or line.startswith("=") or line.endswith(":"):
                continue
            match = LOG_LINE.match(line)
            if not match:
                continue
            label = match.group("label").strip()
            label = re.sub(r"\s+interactions?$", "", label, flags=re.I).strip()
            partners = [p.strip() for p in match.group("partners").split(" and ")]
            # a water bridge names the bridging water in the first partner
            water = ""
            if partners and "," in partners[0]:
                head, _, tail = partners[0].partition(",")
                partners[0] = head.strip()
                water = tail.strip()
            row = {"interaction": label,
                   "partner 1": partners[0] if partners else "",
                   "partner 2": partners[1] if len(partners) > 1 else "",
                   "water": water}
            extra = []
            for piece in match.group("fields").split(","):
                if "=" not in piece:
                    continue
                key, _, value = piece.partition("=")
                key, value = key.strip(), value.strip()
                try:
                    number = round(float(value), 2)
                except ValueError:
                    extra.append("%s = %s" % (key, value))
                    continue
                if "distance" in key.lower():
                    row.setdefault("distance", number)
                elif "angle" in key.lower():
                    row.setdefault("angle", number)
                else:
                    extra.append("%s = %s" % (key, number))
            atoms = match.group("atoms")
            if atoms:
                row["atoms"] = ", ".join(
                    a.strip().strip("'\"") for a in atoms.strip("()").split(",")
                    if a.strip())
            row["other"] = "; ".join(extra)
            rows.append(row)
        columns = ["interaction", "partner 1", "partner 2", "water", "atoms",
                   "distance", "angle", "other"]
        frame = pd.DataFrame(rows)
        for col in columns:
            if col not in frame.columns:
                frame[col] = None
        frame = frame[columns]
        return frame.loc[:, [c for c in columns
                             if frame[c].notna().any() and (frame[c] != "").any()]]

    def to_csv(self, path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        self.table.to_csv(path, index=False)
        return path


def _as_pairs(complex_pdb, mol2, scratch):
    """Resolve the PDB and MOL2 of one complex inside the scratch directory."""
    complex_pdb = Path(complex_pdb)
    if not complex_pdb.is_file():
        raise FileNotFoundError(complex_pdb)
    work = Path(scratch) / "input"
    work.mkdir(parents=True, exist_ok=True)
    pdb_local = work / complex_pdb.name
    shutil.copyfile(complex_pdb, pdb_local)
    if mol2 is None:
        mol2_local = pdb_to_mol2(pdb_local, work / (complex_pdb.stem + ".mol2"))
    else:
        mol2 = Path(mol2)
        if not mol2.is_file():
            raise FileNotFoundError(mol2)
        mol2_local = work / (complex_pdb.stem + ".mol2")
        shutil.copyfile(mol2, mol2_local)
    return work, pdb_local.name, mol2_local.name


def analyse(complex_pdb, mol2=None, save_dir=None, scratch=None, verbose=False,
            capture_log=True):
    """Detect every interaction type in one complex.

    complex_pdb  a docked complex, protein and ligand in one file, with hydrogens
    mol2         the same complex in MOL2, built with Open Babel if left out
    save_dir     write the original per-ligand log and CSV here as well
    scratch      where intermediates go, a fresh temporary directory by default
    capture_log  read the report back onto each result, see InteractionResult.log

    Returns a list of InteractionResult, one for each ligand found in the file.
    """
    owned = scratch is None
    scratch = Path(tempfile.mkdtemp(prefix="poseprint_")) if owned else Path(scratch)
    try:
        work, pdb_name, mol2_name = _as_pairs(complex_pdb, mol2, scratch)
        reports = Path(save_dir) if save_dir else (scratch / "report")
        if save_dir or capture_log:
            reports.mkdir(parents=True, exist_ok=True)
        rows = analysis.run_analyze([pdb_name], [mol2_name], str(work),
                                    save_dir=str(reports) if (save_dir or capture_log)
                                    else None,
                                    scratch=str(scratch / "tmp"), verbose=verbose)
        out = []
        for r in rows:
            text = ""
            if capture_log:
                stem = r["pose"].split(".pdb")[0]
                report = reports / ("results_%s_%s.log" % (stem, r["ligand"]))
                if report.is_file():
                    text = report.read_text(encoding="utf-8", errors="replace")
            out.append(InteractionResult(r["pose"], r["ligand"], r["table"], text))
        return out
    finally:
        if owned:
            shutil.rmtree(scratch, ignore_errors=True)


def analyse_many(complex_pdbs, save_dir=None, verbose=False):
    """Run analyse over several complexes and return the results in one flat list."""
    out = []
    for path in complex_pdbs:
        out.extend(analyse(path, save_dir=save_dir, verbose=verbose))
    return out


def fingerprint(complex_pdb, mol2=None, scratch=None, verbose=False):
    """The hydrogen-bond and hydrophobic summary of one complex.

    This is the reduction the CLINCH labels are built from. One dictionary per ligand,
    carrying the two counts and the residues that make them.
    """
    out = []
    for res in analyse(complex_pdb, mol2=mol2, scratch=scratch, verbose=verbose):
        counts = res.counts()
        out.append({
            "pose": res.pose,
            "ligand": res.ligand,
            "hbond_count": counts.get("H-bonds", 0),
            "hbond_residues": res.residues("H-bonds"),
            "hydrophobic_count": counts.get("Hydrophobic interaction", 0),
            "hydrophobic_residues": res.residues("Hydrophobic interaction"),
        })
    return out
