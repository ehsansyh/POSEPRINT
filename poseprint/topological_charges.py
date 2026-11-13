"""Formal charge assignment by local topology.

Carried over from the working copy of POSEPRINT with its detection logic unchanged.
See build_poseprint_repo.py in the CLINCH working tree for the exact edits that were
applied, which are limited to module naming, relative imports, dead imports and the
removal of one hard-coded scratch path.
"""

# -*- coding: utf-8 -*-
"""Ligand charged groups from connectivity alone, no bond orders and no OpenBabel.

Why this exists. detect_ligand_cation and detect_ligand_anion read
Chem.MolFromPDBFile(...).GetFormalCharge(), but a PDB carries no charge information, so
they only ever saw the ligands whose SMILES happened to be drawn in a charged form.
Measured over the 3148 ligands in this set, that is 5.4 percent, while 54.7 percent carry
a group that is ionised at pH 7.4. PLIP does not read formal charges either, it matches
functional group patterns, which is why it finds salt bridges where this pipeline finds
none.

The patterns below use only the heavy atom graph, which RDKit builds from proximity when
it reads a PDB. No bond order is trusted, nothing is re-perceived, and OpenBabel is not
involved, so none of the bond order corruption that broke the MOL2 route applies here.
Aromaticity, where a rule needs it, comes from ring planarity, the same geometric test
used for pi-stacking.

Known limits, stated so the measurement can check them:
  clean   carboxylate, sulfonate, phosphate, tetrazole, guanidinium, quaternary nitrogen
  harder  plain aliphatic amines, which are the largest group at 37.2 percent and have to
          be separated from amides, anilines and sulfonamides by their neighbourhood
"""
import numpy as np

PLANARITY_CUT = 0.02
AROMATIC_RING_SIZES = (5, 6)


def _ring_coords(mol, conf, ring):
    return np.array([list(conf.GetAtomPosition(i)) for i in ring], dtype=float)


def aromatic_ring_atoms(mol):
    """Indices of atoms sitting in a flat 5 or 6 membered ring."""
    try:
        conf = mol.GetConformer()
    except ValueError:
        return set()
    out = set()
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) not in AROMATIC_RING_SIZES:
            continue
        pts = _ring_coords(mol, conf, ring)
        centred = pts - pts.mean(axis=0)
        normal = np.linalg.svd(centred)[2][-1]
        if float(np.sqrt(np.mean((centred @ normal) ** 2))) <= PLANARITY_CUT:
            out.update(ring)
    return out


def _terminal(atom, symbol):
    return [n for n in atom.GetNeighbors()
            if n.GetSymbol() == symbol and n.GetDegree() == 1]


def anionic_atoms(mol):
    """Atoms carrying a negative charge at pH 7.4, by connectivity pattern."""
    found = {}
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        term_o = _terminal(atom, "O")
        # carboxylate and carboxylic acid. An ester or an anhydride has only one
        # terminal oxygen on each carbon, so neither matches.
        if sym == "C" and len(term_o) >= 2:
            for o in term_o:
                found[o.GetIdx()] = o
        # sulfonate. A sulfonamide carries only two terminal oxygens and is not ionised.
        elif sym == "S" and len(term_o) >= 3:
            for o in term_o:
                found[o.GetIdx()] = o
        # phosphate and phosphonate, including their esters, which stay anionic
        elif sym == "P" and len(term_o) >= 2:
            for o in term_o:
                found[o.GetIdx()] = o

    # tetrazole, a five membered ring of four nitrogens and one carbon
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) != 5:
            continue
        symbols = [mol.GetAtomWithIdx(i).GetSymbol() for i in ring]
        if symbols.count("N") == 4 and symbols.count("C") == 1:
            for i in ring:
                a = mol.GetAtomWithIdx(i)
                if a.GetSymbol() == "N":
                    found[i] = a
    return list(found.values())


def cationic_atoms(mol):
    """Atoms carrying a positive charge at pH 7.4, by connectivity pattern."""
    aromatic = aromatic_ring_atoms(mol)
    found = {}

    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "N":
            continue
        # quaternary nitrogen, permanently charged
        if atom.GetDegree() >= 4:
            found[atom.GetIdx()] = atom

    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "C" or atom.GetIdx() in aromatic:
            continue
        neighbours = list(atom.GetNeighbors())
        nitrogens = [n for n in neighbours if n.GetSymbol() == "N"]
        oxygens = [n for n in neighbours if n.GetSymbol() == "O"]
        if oxygens or len(neighbours) != 3:
            continue
        # three nitrogens is guanidinium, two is amidinium. Both are strong bases.
        # A urea carbon is excluded by the oxygen test, an aromatic ring carbon by the
        # aromatic test, so 2-aminopyrimidine and triazine cores do not match.
        if len(nitrogens) in (2, 3):
            for n in nitrogens:
                found[n.GetIdx()] = n

    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "N" or atom.GetIdx() in aromatic:
            continue
        if atom.GetDegree() > 3:
            continue
        basic = True
        for nb in atom.GetNeighbors():
            if nb.GetSymbol() != "C":
                basic = False           # N-N and N-O are not simple basic amines
                break
            if nb.GetIdx() in aromatic:
                basic = False           # aniline, pKa far below 7.4
                break
            if _terminal(nb, "O"):
                basic = False           # amide
                break
        if basic:
            # a sulfonamide nitrogen is bonded to a sulfur, already excluded above
            found[atom.GetIdx()] = atom
    return list(found.values())


def _named(atoms):
    """The shape the existing detectors return, (index, atom name)."""
    out = []
    for a in atoms:
        info = a.GetPDBResidueInfo()
        name = info.GetName() if info is not None else a.GetSymbol()
        out.append((a.GetIdx(), name))
    return out


def cations(mol):
    return _named(cationic_atoms(mol))


def anions(mol):
    return _named(anionic_atoms(mol))
