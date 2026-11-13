# -*- coding: utf-8 -*-
"""POSEPRINT in the browser.

    streamlit run app/streamlit_app.py

Drop in a docked complex and the page reports which residues contact the ligand and how,
draws the complex with the contacts coloured, and offers the table as a CSV. A second tab
takes a folder's worth of poses and returns one row per pose, which is the form a
scoring or learning model usually wants.
"""
from __future__ import annotations

import sys
import tempfile
import warnings
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
warnings.filterwarnings("ignore")

from poseprint import INTERACTION_COLUMNS, __version__, analyse, have_openbabel  # noqa: E402
from poseprint.diagram import draw_diagram, interaction_records  # noqa: E402

HB = "#08519C"
HY = "#6BAED6"
OTHER = "#D55E00"
PROTEIN = "#D5DBE2"
LIGAND = "#3A3F45"

TYPE_COLOUR = {
    "H-bonds": HB,
    "Hydrophobic interaction": HY,
    "Salt_bridge": OTHER,
    "Water_bridge": OTHER,
    "Pi-Stacking (T-shape)": OTHER,
    "Pi-Stacking(F-F)": OTHER,
    "Pi-Stacking(parallel)": OTHER,
    "Pi-cation": OTHER,
    "Halogen bonds": OTHER,
}

st.set_page_config(page_title="POSEPRINT", page_icon="⬢", layout="wide")


def _store(upload) -> Path:
    tmp = Path(tempfile.mkdtemp(prefix="poseprint_ui_"))
    target = tmp / upload.name
    target.write_bytes(upload.getbuffer())
    return target


@st.cache_data(show_spinner=False)
def _diagram(pdb_bytes: bytes, name: str, title: str, show_distances: bool):
    """Render the diagram once per complex and hand back the three file formats."""
    import io as _io

    import matplotlib.pyplot as plt
    tmp = Path(tempfile.mkdtemp(prefix="poseprint_dia_"))
    pdb = tmp / name
    pdb.write_bytes(pdb_bytes)
    fig = draw_diagram(pdb, title=title or None, show_distances=show_distances)
    files = {}
    for fmt in ("png", "pdf", "svg"):
        buffer = _io.BytesIO()
        fig.savefig(buffer, format=fmt, dpi=600, facecolor="white",
                    bbox_inches="tight", pad_inches=0.04)
        files[fmt] = buffer.getvalue()
    plt.close(fig)
    return files


@st.cache_data(show_spinner=False)
def _run(pdb_bytes: bytes, name: str, mol2_bytes=None, mol2_name=None):
    """Cached so moving a slider does not re-run the detection."""
    tmp = Path(tempfile.mkdtemp(prefix="poseprint_run_"))
    pdb = tmp / name
    pdb.write_bytes(pdb_bytes)
    mol2 = None
    if mol2_bytes is not None:
        mol2 = tmp / (mol2_name or "ligand.mol2")
        mol2.write_bytes(mol2_bytes)
    out = []
    for res in analyse(pdb, mol2=mol2, verbose=False):
        out.append({"ligand": res.ligand, "counts": res.counts(),
                    "table": res.table, "contacts": res.contact_table(),
                    "residues": {c: res.residues(c) for c in INTERACTION_COLUMNS}})
    return out, pdb.read_text(encoding="utf-8", errors="replace")


def _viewer(pdb_text: str, residues: dict, height: int = 520, style: str = "cartoon"):
    import py3Dmol
    view = py3Dmol.view(width="100%", height=height)
    view.addModel(pdb_text, "pdb")
    base = {"cartoon": {"color": PROTEIN, "thickness": 0.6}} if style == "cartoon" \
        else {"line": {"color": PROTEIN}}
    view.setStyle({"hetflag": False}, base)
    view.setStyle({"hetflag": True}, {"stick": {"colorscheme": "default", "radius": 0.18}})
    for column, names in residues.items():
        if not names:
            continue
        colour = TYPE_COLOUR.get(column, OTHER)
        numbers = []
        for label in names:
            parts = str(label).split()
            if len(parts) == 2 and parts[1].lstrip("-").isdigit():
                numbers.append(int(parts[1]))
        if numbers:
            view.addStyle({"resi": numbers, "hetflag": False},
                          {"stick": {"color": colour, "radius": 0.16}})
            view.addStyle({"resi": numbers, "hetflag": False},
                          {"cartoon": {"color": colour}})
    # fit the whole complex, so the pocket is seen in context. The viewer is
    # interactive, so anyone who wants the pocket can scroll into it.
    view.zoomTo()
    return view


def panel_single():
    left, right = st.columns([1, 1.25])

    with left:
        st.subheader("Complex")
        example = Path(__file__).resolve().parent.parent / "examples" / \
            "hiv1_protease_darunavir.pdb"
        use_example = st.checkbox("Use the bundled example, HIV-1 protease with darunavir",
                                  value=example.is_file())
        upload = None
        mol2_up = None
        if not use_example:
            upload = st.file_uploader("A docked complex in PDB format, with hydrogens",
                                      type=["pdb"])
            with st.expander("Supply the MOL2 yourself, optional"):
                st.caption("Left out, the MOL2 is built with Open Babel. It is only used "
                           "for ring perception, so the hydrogen-bond and hydrophobic "
                           "results do not depend on it.")
                mol2_up = st.file_uploader("MOL2 of the same complex", type=["mol2"])

        if use_example and example.is_file():
            pdb_bytes, pdb_name = example.read_bytes(), example.name
        elif upload is not None:
            pdb_bytes, pdb_name = upload.getbuffer().tobytes(), upload.name
        else:
            st.info("Choose a complex to begin.")
            return

        if not have_openbabel() and mol2_up is None:
            st.warning("Open Babel is not on the path. The MOL2 will be built with "
                       "MDAnalysis, which guesses bonds, so the aromatic interactions "
                       "may differ. Hydrogen bonds and hydrophobic contacts are "
                       "unaffected.")

        style = st.radio("Protein style", ["cartoon", "lines"], horizontal=True)

    with st.spinner("Detecting interactions"):
        try:
            results, pdb_text = _run(bytes(pdb_bytes), pdb_name,
                                     mol2_up.getbuffer().tobytes() if mol2_up else None,
                                     mol2_up.name if mol2_up else None)
        except Exception as exc:
            st.error("POSEPRINT could not read this complex. %s: %s"
                     % (type(exc).__name__, exc))
            return

    if not results:
        st.warning("No ligand was found. POSEPRINT reads the ligand from the HETATM "
                   "records, so the complex needs the ligand in the same file.")
        return

    names = [r["ligand"] for r in results]
    chosen = names[0] if len(names) == 1 else st.selectbox("Ligand", names)
    res = next(r for r in results if r["ligand"] == chosen)

    with right:
        st.subheader("Structure")
        import streamlit.components.v1 as components
        view = _viewer(pdb_text, res["residues"], style=style)
        components.html(view._make_html(), height=540, scrolling=False)
        st.caption("Ligand by element. Hydrogen-bond partners in dark blue, hydrophobic "
                   "contacts in light blue, every other interaction type in orange.")

    counts = {k: v for k, v in res["counts"].items() if v}
    if not counts:
        st.info("No interaction of any type was detected for %s. If the complex has no "
                "hydrogens, that is the usual reason." % chosen)
    else:
        cols = st.columns(min(len(counts), 5))
        for (name, value), col in zip(counts.items(), cols * 5):
            col.metric(name.replace("_", " "), value)

    st.subheader("Interaction diagram")
    st.caption("The ligand with the residues it touches. Dashed lines run to the "
               "ligand atom the interaction is made with. A spoked arc is a hydrophobic "
               "contact, which the detector reports per residue rather than per atom, "
               "so it faces the nearest ligand atom. Residues are coloured by class.")
    opt_left, opt_right = st.columns([2, 1])
    with opt_left:
        diagram_title = st.text_input("Caption drawn on the figure, leave empty for none",
                                      value=Path(pdb_name).stem)
    with opt_right:
        with_distances = st.checkbox("Mark the distances", value=True)
    try:
        files = _diagram(bytes(pdb_bytes), pdb_name, diagram_title, with_distances)
    except Exception as exc:
        st.error("The diagram could not be drawn. %s: %s" % (type(exc).__name__, exc))
    else:
        st.image(files["png"], use_container_width=True)
        stem = "poseprint_%s_%s" % (Path(pdb_name).stem, chosen)
        d1, d2, d3 = st.columns(3)
        d1.download_button("PDF for the manuscript", files["pdf"],
                           file_name=stem + ".pdf", mime="application/pdf",
                           use_container_width=True)
        d2.download_button("SVG to edit", files["svg"], file_name=stem + ".svg",
                           mime="image/svg+xml", use_container_width=True)
        d3.download_button("PNG", files["png"], file_name=stem + ".png",
                           mime="image/png", use_container_width=True)
        st.caption("The PDF and the SVG are vector, so they stay sharp at any size in a "
                   "manuscript. The SVG opens in Illustrator or Inkscape if a label "
                   "needs moving.")

    st.subheader("Residues in contact")
    table = res["contacts"]
    st.dataframe(table, use_container_width=True, hide_index=True)
    st.download_button("Download the full table as CSV",
                       res["table"].to_csv(index=False).encode("utf-8"),
                       file_name="poseprint_%s_%s.csv" % (Path(pdb_name).stem, chosen),
                       mime="text/csv")


def panel_batch():
    st.subheader("Many poses at once")
    st.caption("One row per pose and ligand, with the counts and the residues that make "
               "them. This is the form the CLINCH labels were built in.")
    uploads = st.file_uploader("Complexes in PDB format", type=["pdb"],
                               accept_multiple_files=True)
    if not uploads:
        st.info("Add one or more files to begin.")
        return
    if st.button("Run %d complex(es)" % len(uploads), type="primary"):
        rows, failed = [], []
        bar = st.progress(0.0, text="Starting")
        for i, up in enumerate(uploads, 1):
            bar.progress(i / len(uploads), text="%s (%d of %d)" % (up.name, i, len(uploads)))
            try:
                results, _ = _run(up.getbuffer().tobytes(), up.name)
            except Exception as exc:
                failed.append((up.name, "%s: %s" % (type(exc).__name__, exc)))
                continue
            for res in results:
                rows.append({
                    "complex": up.name,
                    "ligand": res["ligand"],
                    "hbond_count": res["counts"].get("H-bonds", 0),
                    "hbond_residues": ";".join(res["residues"]["H-bonds"]),
                    "hydrophobic_count": res["counts"].get("Hydrophobic interaction", 0),
                    "hydrophobic_residues": ";".join(
                        res["residues"]["Hydrophobic interaction"]),
                    **{c: res["counts"].get(c, 0) for c in INTERACTION_COLUMNS
                       if c not in ("H-bonds", "Hydrophobic interaction")},
                })
        bar.empty()
        if rows:
            frame = pd.DataFrame(rows)
            st.dataframe(frame, use_container_width=True, hide_index=True)
            st.download_button("Download as CSV",
                               frame.to_csv(index=False).encode("utf-8"),
                               file_name="poseprint_batch.csv", mime="text/csv")
        if failed:
            st.warning("%d complex(es) could not be read." % len(failed))
            for name, why in failed:
                st.text("  %s  %s" % (name, why))


st.title("POSEPRINT")
st.caption("Protein-ligand interaction fingerprints from docked poses, version %s"
           % __version__)

tab_one, tab_many, tab_about = st.tabs(["One complex", "A batch", "About"])
with tab_one:
    panel_single()
with tab_many:
    panel_batch()
with tab_about:
    st.markdown(
        """
POSEPRINT reads a docked complex and reports, residue by residue, which of seven
interaction types it makes with the ligand.

**The seven types**, hydrogen bonds, hydrophobic contacts, salt bridges, water bridges,
pi stacking in its parallel, offset and T-shaped forms, pi-cation interactions and
halogen bonds.

**What the complex needs.** Protein and ligand in one PDB file, the ligand written as
HETATM records, and hydrogens present. Hydrogen bonds are found through the donor
hydrogen, so a heavy-atom-only file returns nothing rather than failing.

**How the counts are read.** The table carries one row per residue and one row for the
ligand. Every contact appears twice, once against the residue and once against the
ligand, so the ligand row is the total. A residue label carries no chain, so in a
homodimer the two copies of a residue share one label.

**Citation.** POSEPRINT was written for the CLINCH study of HIV-1 protease inhibitor
design, where it labelled 142,160 docked poses of 3,949 ligands. Nothing in it is
specific to that target.
        """)
