# -*- coding: utf-8 -*-
"""Build the example complex that ships with the repository.

The source is the prepared wild-type HIV-1 protease and darunavir complex that entered
the molecular dynamics in the CLINCH study. Two things have to be right for a complex to
be readable by POSEPRINT.

  Hydrogens have to be present. Hydrogen bonds are found through the donor hydrogen, so a
  heavy-atom-only file silently returns nothing.

  Every atom in a residue needs its own name. Maestro does not always carry PDB atom
  names for a ligand, and a PDB file where seventy-five ligand atoms are all called C or
  H collapses to five atoms when any PDB parser reads it, because a parser keeps one atom
  per name per residue. Names are generated here where the source leaves them blank.

  python make_example.py
"""
from __future__ import annotations

import collections
import io
import shlex
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
HERE = Path(__file__).resolve().parent
SRC = Path(r"C:/Users/ehsan/Desktop/phd-codes/CLINCH_H_MD_200ns/HIV_WT_DARUNAVIR.mae")
DST = HERE / "hiv1_protease_darunavir.pdb"

SYM = {1: "H", 6: "C", 7: "N", 8: "O", 9: "F", 11: "NA", 12: "MG", 15: "P",
       16: "S", 17: "CL", 19: "K", 20: "CA", 30: "ZN", 35: "BR", 53: "I"}


def convert(src: Path, dst: Path) -> dict:
    txt = src.read_text(encoding="utf-8", errors="replace")
    i = txt.index("m_atom[")
    head_end = txt.index(":::", i)
    props = [l.strip() for l in txt[txt.index("{", i) + 1:head_end].split("\n")
             if l.strip() and not l.strip().startswith("#")]
    body = txt[head_end + 3:txt.index(":::", head_end + 3)]
    tok = shlex.split(body)
    ncol = len(props) + 1
    idx = {name: props.index(name) + 1 for name in props}

    lines, serial = [], 0
    used = collections.defaultdict(set)      # per residue, the names already taken
    counter = collections.defaultdict(int)   # per residue and element, a running index
    chains, residues, renamed = set(), set(), 0

    for k in range(0, len(tok) - ncol + 1, ncol):
        row = tok[k:k + ncol]
        z = int(row[idx["i_m_atomic_number"]])
        el = SYM.get(z, "C")
        resn = (row[idx["s_m_pdb_residue_name"]] or "UNK").strip()[:3] or "UNK"
        chain = (row[idx["s_m_chain_name"]] or "A").strip()[:1] or "A"
        resi = int(row[idx["i_m_residue_number"]])
        name = (row[idx["s_m_pdb_atom_name"]] or "").strip()[:4]
        key = (chain, resi, resn)

        if not name or name in used[key]:
            while True:
                counter[(key, el)] += 1
                cand = ("%s%d" % (el, counter[(key, el)]))[:4]
                if cand not in used[key]:
                    break
            if name:
                renamed += 1
            name = cand
        used[key].add(name)

        serial += 1
        rec = "HETATM" if resn == "UNK" else "ATOM  "
        lines.append("%-6s%5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s"
                     % (rec, serial % 100000, name, resn, chain, resi,
                        float(row[idx["r_m_x_coord"]]), float(row[idx["r_m_y_coord"]]),
                        float(row[idx["r_m_z_coord"]]), 1.00, 0.00, el.rjust(2)))
        chains.add(chain)
        residues.add(key)

    lines.append("END")
    dst.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return {"atoms": serial, "chains": sorted(chains), "residues": len(residues),
            "renamed": renamed}


def main() -> None:
    if not SRC.is_file():
        print("  the source complex is not on this machine, %s" % SRC)
        sys.exit(1)
    info = convert(SRC, DST)
    print("  %s -> %s" % (SRC.name, DST.name))
    print("  %d atoms, chains %s, %d residues, %d atom names generated"
          % (info["atoms"], info["chains"], info["residues"], info["renamed"]))

    from Bio.PDB import PDBParser
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        st = PDBParser(QUIET=True).get_structure("c", str(DST))
    kept = len(list(st.get_atoms()))
    lig = [r for r in st.get_residues() if r.get_id()[0] != " "]
    print("  a PDB parser keeps %d of %d atoms, ligand residues %s"
          % (kept, info["atoms"],
             [(r.get_resname(), len(list(r.get_atoms()))) for r in lig]))
    if kept != info["atoms"]:
        print("  WARNING atoms are still being dropped, names are not unique")
        sys.exit(1)


if __name__ == "__main__":
    main()
