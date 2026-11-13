# -*- coding: utf-8 -*-
"""A publication-ready two-dimensional interaction diagram.

The ligand is drawn as a structure, the residues that touch it sit around it as coloured
nodes, and every interaction that names a ligand atom is drawn as a dashed line to that
atom. Hydrophobic contacts carry no atom in the detector output, so they are drawn the
way LigPlot draws them, as a spoked arc facing the ligand.

    from poseprint import draw_diagram
    draw_diagram("complex.pdb", out="figure.pdf")

Residues are coloured by class, which is the thing a reader wants to see first.

    nonpolar   Ala Val Leu Ile Pro Met Phe Trp Gly
    polar      Ser Thr Cys Asn Gln Tyr
    acidic     Asp Glu
    basic      Lys Arg His
"""
from __future__ import annotations

import io
import math
import shutil
import tempfile
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
from PIL import Image

__all__ = ["draw_diagram", "interaction_records", "RESIDUE_CLASS", "CLASS_COLOUR"]

# --------------------------------------------------------------------- styling
RESIDUE_CLASS = {
    **{r: "nonpolar" for r in
       ("ALA", "VAL", "LEU", "ILE", "PRO", "MET", "PHE", "TRP", "GLY")},
    **{r: "polar" for r in ("SER", "THR", "CYS", "ASN", "GLN", "TYR")},
    **{r: "acidic" for r in ("ASP", "GLU")},
    **{r: "basic" for r in ("LYS", "ARG", "HIS")},
}

CLASS_COLOUR = {
    "nonpolar": ("#D6E8C8", "#4E7A3A"),
    "polar": ("#CFE0EF", "#2E6D9E"),
    "acidic": ("#F6CFC9", "#B03A2E"),
    "basic": ("#DDD2EC", "#6A4C93"),
    "other": ("#E6E9EC", "#79838C"),
}

INTERACTION_STYLE = {
    "hydrogen bond": ("#08519C", (0, (3.0, 1.8))),
    "water bridge": ("#4FA3D1", (0, (1.4, 1.4))),
    "salt bridge": ("#D55E00", (0, (4.5, 1.8))),
    "halogen bond": ("#0E8F7E", (0, (3.0, 1.2, 1.0, 1.2))),
    "pi stacking": ("#8E44AD", (0, (2.4, 1.6))),
    "pi-cation": ("#C2278F", (0, (2.4, 1.6))),
}
HYDROPHOBIC_COLOUR = "#7E9E5E"

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "mathtext.fontset": "stix",
})
TIMES = Path("C:/Windows/Fonts/times.ttf")


def _title(resname: str) -> str:
    """ASP -> Asp, the way a residue is written in a paper."""
    return resname[:1].upper() + resname[1:].lower()


def residue_class(resname: str) -> str:
    return RESIDUE_CLASS.get(resname.upper(), "other")


# ------------------------------------------------------------------- detection
def interaction_records(complex_pdb, mol2=None, scratch=None):
    """Every interaction, one record each, with the ligand atom where it is known.

    Returns a dict with the ligand identity and a list of records of the form
    {type, resname, resid, ligand_atom, distance}. The detectors report an atom pair for
    hydrogen bonds, water bridges and halogen bonds. Hydrophobic contacts, salt bridges
    and the aromatic interactions are reported per residue only, so ligand_atom is None
    and the diagram anchors them to the ligand as a whole.
    """
    from Bio.PDB import PDBParser

    from . import halogen, hbond, hydrophobic, pi_cation, pi_stacking
    from . import salt_bridge as sb
    from . import water_bridge as wb
    from .analysis import ext_lig
    from .convert import pdb_to_mol2

    owned = scratch is None
    scratch = Path(tempfile.mkdtemp(prefix="poseprint_dia_")) if owned else Path(scratch)
    try:
        work = scratch / "in"
        work.mkdir(parents=True, exist_ok=True)
        pdb_local = work / Path(complex_pdb).name
        shutil.copyfile(complex_pdb, pdb_local)
        mol2_local = (work / (pdb_local.stem + ".mol2"))
        if mol2 is None:
            pdb_to_mol2(pdb_local, mol2_local)
        else:
            shutil.copyfile(mol2, mol2_local)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            structure = PDBParser(QUIET=True).get_structure("c", str(pdb_local))
        ligands, ligands_mol2 = ext_lig(structure)
        if not ligands:
            return []
        dirct = str(scratch / "tmp") + "/"
        Path(dirct).mkdir(parents=True, exist_ok=True)

        out = []
        for lig_id, lig_id_mol2 in zip(ligands, ligands_mol2):
            lig_resname = lig_id_mol2.split()[0]
            records = []

            def _atom_pairs(table, kind, atom_index):
                """The detectors that key on a residue pair and carry atom names."""
                for key, values in (table or {}).items():
                    if len(key) < 4:
                        continue
                    a_name, a_id, b_name, b_id = key[0], key[1], key[2], key[3]
                    lig_first = a_name == lig_resname
                    if not lig_first and b_name != lig_resname:
                        continue
                    res_name, res_id = (b_name, b_id) if lig_first else (a_name, a_id)
                    if res_name == lig_resname:
                        continue
                    for entry in values:
                        atoms = entry[0] if isinstance(entry, (tuple, list)) else None
                        dist = entry[1] if len(entry) > 1 else None
                        lig_atom = None
                        if atoms and len(atoms) > max(atom_index):
                            lig_atom = atoms[atom_index[0] if lig_first else atom_index[1]]
                        records.append({"type": kind, "resname": res_name,
                                        "resid": int(res_id), "ligand_atom": lig_atom,
                                        "distance": float(dist) if dist else None})

            try:
                _atom_pairs(hbond.main(str(pdb_local), structure), "hydrogen bond", (0, 2))
            except Exception:
                pass
            try:
                _atom_pairs(wb.main(str(pdb_local), structure), "water bridge", (0, 2))
            except Exception:
                pass
            try:
                _atom_pairs(halogen.main(str(pdb_local), structure), "halogen bond", (0, 1))
            except Exception:
                pass
            try:
                for row in hydrophobic.hydrophobic(str(pdb_local), str(mol2_local),
                                                   lig_id, lig_id_mol2, dirct,
                                                   structure) or []:
                    if len(row) >= 4:
                        records.append({"type": "hydrophobic", "resname": row[2],
                                        "resid": int(row[3]), "ligand_atom": None,
                                        "distance": float(row[4]) if len(row) > 4 else None})
            except Exception:
                pass
            try:
                for pair in sb.run(str(pdb_local), str(mol2_local), lig_id,
                                   lig_id_mol2, dirct, structure) or []:
                    for atom in pair:
                        parent = atom.get_parent()
                        if parent.get_resname() != lig_resname:
                            records.append({"type": "salt bridge",
                                            "resname": parent.get_resname(),
                                            "resid": int(parent.get_id()[1]),
                                            "ligand_atom": None, "distance": None})
                            break
            except Exception:
                pass
            try:
                tshape, parallel, ff = pi_stacking.run(str(mol2_local), str(pdb_local),
                                                       lig_id, lig_id_mol2, dirct,
                                                       structure)
                for group in (tshape, parallel, ff):
                    for key in group or []:
                        for item in key:
                            if hasattr(item, "get_parent"):
                                parent = item.get_parent()
                                if parent.get_resname() != lig_resname:
                                    records.append({"type": "pi stacking",
                                                    "resname": parent.get_resname(),
                                                    "resid": int(parent.get_id()[1]),
                                                    "ligand_atom": None,
                                                    "distance": None})
                                    break
            except Exception:
                pass
            try:
                lig_cat, prot_cat = pi_cation.run(str(pdb_local), str(mol2_local), lig_id,
                                                  lig_id_mol2, dirct, structure)
                for group in (lig_cat, prot_cat):
                    for key in group or []:
                        for item in key:
                            if hasattr(item, "get_parent"):
                                parent = item.get_parent()
                                if parent.get_resname() != lig_resname:
                                    records.append({"type": "pi-cation",
                                                    "resname": parent.get_resname(),
                                                    "resid": int(parent.get_id()[1]),
                                                    "ligand_atom": None,
                                                    "distance": None})
                                    break
            except Exception:
                pass

            # a hydrophobic contact or a salt bridge is reported per residue, so the
            # diagram would have nothing to point at. The nearest ligand heavy atom to
            # that residue is a measured fact, not a guess, and it places the node
            # where the contact actually is.
            need = {(r["resname"], r["resid"]) for r in records if not r["ligand_atom"]}
            if need:
                lig_atoms = [a for res in structure.get_residues()
                             if res.get_resname() == lig_resname
                             for a in res.get_atoms() if a.element != "H"]
                nearest = {}
                for res in structure.get_residues():
                    key = (res.get_resname(), int(res.get_id()[1]))
                    if key not in need or not lig_atoms:
                        continue
                    best, best_d = None, 1e9
                    for atom in res.get_atoms():
                        if atom.element == "H":
                            continue
                        for lat in lig_atoms:
                            d = float(np.linalg.norm(atom.coord - lat.coord))
                            if d < best_d:
                                best, best_d = lat.get_name(), d
                    if best:
                        nearest[key] = best
                for rec in records:
                    if not rec["ligand_atom"]:
                        rec["nearest_atom"] = nearest.get((rec["resname"], rec["resid"]))

            out.append({"ligand": lig_id, "ligand_resname": lig_resname,
                        "pdb": str(complex_pdb), "records": records})
        return out
    finally:
        if owned:
            shutil.rmtree(scratch, ignore_errors=True)


# ---------------------------------------------------------------- ligand image
def _ligand_mol(complex_pdb, lig_resname):
    """The ligand as an RDKit molecule with its PDB atom names and real bond orders."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.*")
    lines = [l for l in Path(complex_pdb).read_text(errors="replace").splitlines(True)
             if l.startswith("HETATM") and l[17:20].strip() == lig_resname]
    if not lines:
        raise ValueError("no HETATM records for %s" % lig_resname)
    mol = Chem.MolFromPDBBlock("".join(lines) + "END\n", removeHs=False, sanitize=False)
    if mol is None:
        raise ValueError("the ligand could not be read from the PDB records")
    try:
        from rdkit.Chem import rdDetermineBonds
        rw = Chem.RWMol(mol)
        rdDetermineBonds.DetermineBondOrders(rw, charge=0)
        mol = rw.GetMol()
    except Exception:
        pass                                   # single bonds still draw a usable sketch
    try:
        mol = Chem.RemoveHs(mol, sanitize=False)
    except Exception:
        pass
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_SETAROMATICITY |
                         Chem.SanitizeFlags.SANITIZE_SYMMRINGS, catchErrors=True)
    AllChem.Compute2DCoords(mol)
    return mol


def _render(mol, px, line_mm, font_mm, box_mm):
    """Draw the ligand and report where every named atom landed, in image pixels."""
    from rdkit.Chem.Draw import rdMolDraw2D
    ppm = px[0] / box_mm
    drawer = rdMolDraw2D.MolDraw2DCairo(px[0], px[1])
    opts = drawer.drawOptions()
    if TIMES.is_file():
        opts.fontFile = str(TIMES)
    opts.bondLineWidth = line_mm * ppm
    opts.fixedFontSize = int(round(font_mm * ppm))
    opts.minFontSize = opts.fixedFontSize
    opts.clearBackground = False
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    image = Image.open(io.BytesIO(drawer.GetDrawingText())).convert("RGBA")
    coords = {}
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None:
            continue
        point = drawer.GetDrawCoords(atom.GetIdx())
        coords[info.GetName().strip()] = (point.x, point.y)
    return image, coords


# -------------------------------------------------------------------- the plot
def _relax(start, sizes, keep_out, bounds, rounds=1200):
    """Settle the labels so none of them overlap each other or the ligand.

    Each label wants to stay where its contact is, so it is pulled back towards its
    starting point every round. The pushes that separate two labels, or push one off
    the structure, are stronger, so collisions always win over the preference.
    """
    pos = np.array(start, dtype=float)
    want = pos.copy()
    half = np.array(sizes, dtype=float) / 2.0
    ox0, oy0, ox1, oy1 = keep_out
    bx0, by0, bx1, by1 = bounds
    for _ in range(rounds):
        shift = np.zeros_like(pos)
        for i in range(len(pos)):
            for j in range(i + 1, len(pos)):
                delta = pos[j] - pos[i]
                need = half[i] + half[j] + np.array([2.4, 1.8])
                overlap = need - np.abs(delta)
                if overlap[0] > 0 and overlap[1] > 0:
                    axis = 0 if overlap[0] < overlap[1] else 1
                    push = overlap[axis] / 2.0 + 0.05
                    sign = 1.0 if delta[axis] >= 0 else -1.0
                    shift[j][axis] += push * sign
                    shift[i][axis] -= push * sign
        for i in range(len(pos)):
            # keep clear of the drawn structure
            x0, y0 = pos[i] - half[i] - np.array([1.6, 1.2])
            x1, y1 = pos[i] + half[i] + np.array([1.6, 1.2])
            if x1 > ox0 and x0 < ox1 and y1 > oy0 and y0 < oy1:
                options = [(ox0 - x1, 0), (ox1 - x0, 0), (0, oy0 - y1), (0, oy1 - y0)]
                dx, dy = min(options, key=lambda d: abs(d[0]) + abs(d[1]))
                shift[i] += np.array([dx, dy])
        pos += shift * 0.75 + (want - pos) * 0.02
        pos[:, 0] = np.clip(pos[:, 0], bx0 + half[:, 0], bx1 - half[:, 0])
        pos[:, 1] = np.clip(pos[:, 1], by0 + half[:, 1], by1 - half[:, 1])
        if np.abs(shift).max() < 0.02:
            break
    return pos


def _hydrophobic_arc(ax, node, towards, colour, clearance, spokes=9):
    """The spoked arc LigPlot uses for a hydrophobic contact.

    It sits just outside the residue label on the side facing the ligand, with the
    spokes pointing at the ligand, which is how the convention reads.
    """
    base = math.atan2(towards[1] - node[1], towards[0] - node[0])
    radius = clearance + 0.9
    span = math.radians(94)
    arc = np.linspace(base - span / 2, base + span / 2, 48)
    ax.plot(node[0] + radius * np.cos(arc), node[1] + radius * np.sin(arc),
            color=colour, lw=1.0, zorder=3, solid_capstyle="round")
    for k in range(spokes):
        angle = base - span / 2 + span * k / (spokes - 1)
        x0 = node[0] + radius * math.cos(angle)
        y0 = node[1] + radius * math.sin(angle)
        x1 = node[0] + (radius + 1.5) * math.cos(angle)
        y1 = node[1] + (radius + 1.5) * math.sin(angle)
        ax.plot([x0, x1], [y0, y1], color=colour, lw=0.8, zorder=3,
                solid_capstyle="round")


def draw_diagram(complex_pdb, out=None, mol2=None, ligand=None, title=None,
                 width_mm=183.0, dpi=600, show_distances=True, records=None):
    """Draw the interaction diagram and optionally save it.

    out may be a path with any extension matplotlib writes, so .pdf and .svg give a
    vector figure that scales into a manuscript without going soft.

    Returns the matplotlib Figure, so a caller can alter it before saving.
    """
    if records is None:
        found = interaction_records(complex_pdb, mol2=mol2)
    else:
        found = records
    if not found:
        raise ValueError("no ligand was found in %s" % complex_pdb)
    entry = found[0] if ligand is None else next(
        (f for f in found if f["ligand"] == ligand), found[0])

    mol = _ligand_mol(complex_pdb, entry["ligand_resname"])

    # one entry per residue, carrying every interaction it makes
    by_residue = {}
    for rec in entry["records"]:
        key = (rec["resname"], rec["resid"])
        slot = by_residue.setdefault(key, {"types": {}, "atoms": set(), "near": set()})
        slot["types"].setdefault(rec["type"], []).append(rec)
        if rec["ligand_atom"]:
            slot["atoms"].add(rec["ligand_atom"])
        elif rec.get("nearest_atom"):
            slot["near"].add(rec["nearest_atom"])

    mm = 1 / 25.4
    height_mm = width_mm * 0.64
    lig_box = width_mm * 0.54
    fig = plt.figure(figsize=(width_mm * mm, height_mm * mm))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, width_mm)
    ax.set_ylim(height_mm, 0)
    ax.axis("off")

    px = (2600, int(2600 * 0.78))
    image, coords = _render(mol, px, line_mm=0.16, font_mm=1.25, box_mm=lig_box)
    array = np.asarray(image)
    ys, xs = np.where(array[..., 3] > 8)
    crop = (xs.min(), ys.min(), xs.max() + 1, ys.max() + 1)
    cropped = image.crop(crop)
    scale = min(lig_box / cropped.size[0], (height_mm * 0.50) / cropped.size[1])
    draw_w, draw_h = cropped.size[0] * scale, cropped.size[1] * scale
    origin = (width_mm / 2 - draw_w / 2, height_mm / 2 - draw_h / 2)
    target = (max(1, round(draw_w / 25.4 * dpi)), max(1, round(draw_h / 25.4 * dpi)))
    inset = fig.add_axes([origin[0] / width_mm, 1 - (origin[1] + draw_h) / height_mm,
                          draw_w / width_mm, draw_h / height_mm])
    inset.imshow(np.asarray(cropped.resize(target, Image.LANCZOS)), interpolation="none")
    inset.axis("off")
    inset.patch.set_alpha(0)

    def to_mm(name):
        if name not in coords:
            return None
        x_px, y_px = coords[name]
        return (origin[0] + (x_px - crop[0]) * scale,
                origin[1] + (y_px - crop[1]) * scale)

    centre = (width_mm / 2, height_mm / 2)
    keys = list(by_residue)
    anchors = []
    for key in keys:
        slot = by_residue[key]
        source = slot["atoms"] or slot["near"]
        points = [p for p in (to_mm(a) for a in source) if p]
        anchors.append((float(np.mean([p[0] for p in points])),
                        float(np.mean([p[1] for p in points]))) if points else centre)

    labels = ["%s%d" % (_title(k[0]), k[1]) for k in keys]
    sizes = [(2.4 + 1.78 * len(t), 6.4) for t in labels]
    structure_box = (origin[0] - 1.0, origin[1] - 1.0,
                     origin[0] + draw_w + 1.0, origin[1] + draw_h + 1.0)

    # Several residues often touch the same end of the ligand, so their contacts give
    # nearly the same direction and they would all start in one spot. The angles are
    # spread around the ring first, keeping their order, and the relaxation below only
    # has to clean up what is left.
    n = len(keys)
    raw = [math.atan2(a[1] - centre[1], a[0] - centre[0]) for a in anchors]
    order = sorted(range(n), key=lambda i: raw[i])
    spread = list(raw)
    gap = 2 * math.pi / max(n, 1) * 0.92
    for _ in range(400):
        moved = False
        for k in range(n):
            i, j = order[k], order[(k + 1) % n]
            delta = (spread[j] - spread[i]) % (2 * math.pi)
            if delta < gap:
                push = (gap - delta) / 2.0
                spread[i] -= push
                spread[j] += push
                moved = True
        if not moved:
            break

    widest = max(w for w, _ in sizes)
    ring_x = max(draw_w / 2 + widest / 2 + 7.0, width_mm * 0.30)
    ring_y = max(draw_h / 2 + 7.0, height_mm * 0.30)
    start = [(centre[0] + ring_x * math.cos(t), centre[1] + ring_y * math.sin(t))
             for t in spread]
    margin = 2.0
    nodes = _relax(start, sizes, structure_box,
                   (margin, margin + 5.0, width_mm - margin, height_mm - margin - 7.0))

    for i, key in enumerate(keys):
        resname, _ = key
        slot = by_residue[key]
        label, (bw, bh) = labels[i], sizes[i]
        fill, edge = CLASS_COLOUR[residue_class(resname)]
        node = (float(nodes[i][0]), float(nodes[i][1]))

        for kind, recs in slot["types"].items():
            if kind == "hydrophobic":
                continue
            colour, dash = INTERACTION_STYLE.get(kind, ("#555555", (0, (3, 2))))
            seen = set()
            for rec in recs:
                target_pt = to_mm(rec["ligand_atom"]) if rec["ligand_atom"] else anchors[i]
                if target_pt is None or target_pt in seen:
                    continue
                seen.add(target_pt)
                ax.plot([node[0], target_pt[0]], [node[1], target_pt[1]], color=colour,
                        lw=0.9, ls=dash, zorder=2, solid_capstyle="round")
                if show_distances and rec["distance"]:
                    t = 0.58
                    mid = (node[0] + (target_pt[0] - node[0]) * t,
                           node[1] + (target_pt[1] - node[1]) * t)
                    ax.text(mid[0], mid[1], "%.2f" % rec["distance"], fontsize=4.8,
                            color=colour, ha="center", va="center", zorder=5,
                            bbox=dict(boxstyle="round,pad=0.14", fc="white", ec="none",
                                      alpha=0.9))
        if "hydrophobic" in slot["types"]:
            _hydrophobic_arc(ax, node, anchors[i], HYDROPHOBIC_COLOUR,
                             clearance=max(bw, bh) / 2.0)

        ax.add_patch(FancyBboxPatch((node[0] - bw / 2, node[1] - bh / 2), bw, bh,
                                    boxstyle="round,pad=0,rounding_size=1.5",
                                    linewidth=0.9, edgecolor=edge, facecolor=fill,
                                    zorder=4, mutation_aspect=1))
        ax.text(node[0], node[1], label, fontsize=6.8, color="#15181C", ha="center",
                va="center", zorder=5)

    if title:
        ax.text(width_mm / 2, 4.0, title, fontsize=8.0, ha="center", va="top",
                color="#15181C")

    used_types = sorted({k for s in by_residue.values() for k in s["types"]})
    handles = []
    for kind in used_types:
        if kind == "hydrophobic":
            handles.append(Line2D([0], [0], color=HYDROPHOBIC_COLOUR, lw=1.0,
                                  label="hydrophobic contact"))
        else:
            colour, dash = INTERACTION_STYLE.get(kind, ("#555555", (0, (3, 2))))
            handles.append(Line2D([0], [0], color=colour, lw=1.0, ls=dash, label=kind))
    classes = sorted({residue_class(k[0]) for k in by_residue})
    for klass in classes:
        fill, edge = CLASS_COLOUR[klass]
        handles.append(Line2D([0], [0], marker="s", lw=0, markersize=5.5,
                              markerfacecolor=fill, markeredgecolor=edge,
                              label=klass))
    if handles:
        ax.legend(handles=handles, loc="lower center", ncol=min(len(handles), 6),
                  frameon=False, fontsize=5.8, handlelength=2.2,
                  bbox_to_anchor=(0.5, -0.012), columnspacing=1.4)

    if out:
        out = Path(out)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=dpi, facecolor="white", bbox_inches="tight",
                    pad_inches=0.04)
    return fig
