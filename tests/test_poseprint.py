# -*- coding: utf-8 -*-
"""Tests that pin the chemistry, not just the plumbing.

The example is darunavir in wild-type HIV-1 protease. Its binding mode is settled, so
the contacts below are a real regression test. If a change to a detector moves them, the
change is wrong until proven otherwise.
"""
from __future__ import annotations

import warnings
from pathlib import Path

import pytest

from poseprint import analyse, fingerprint
from poseprint.core import INTERACTION_COLUMNS

warnings.filterwarnings("ignore")
EXAMPLE = Path(__file__).resolve().parent.parent / "examples" / \
    "hiv1_protease_darunavir.pdb"
pytestmark = pytest.mark.skipif(not EXAMPLE.is_file(),
                                reason="the example complex is not in the checkout")


@pytest.fixture(scope="module")
def result():
    found = analyse(EXAMPLE)
    assert len(found) == 1, "the example carries exactly one ligand"
    return found[0]


def test_ligand_is_found(result):
    assert result.ligand.startswith("UNK")


def test_catalytic_aspartate_is_a_hydrogen_bond_partner(result):
    """Asp25 is the catalytic residue every HIV-1 protease inhibitor engages."""
    partners = result.residues("H-bonds")
    assert "ASP 25" in partners, partners


def test_hydrogen_bond_partners_are_exactly_the_known_pair(result):
    assert sorted(result.residues("H-bonds")) == ["ASP 25", "ASP 30"]


def test_hydrophobic_contacts_are_the_canonical_pocket(result):
    """The S1, S1', S2 and S2' subsites darunavir is known to occupy."""
    assert sorted(result.residues("Hydrophobic interaction")) == [
        "ALA 28", "ILE 47", "ILE 54", "ILE 84", "LEU 23", "VAL 82"]


def test_counts_come_from_the_ligand_row(result):
    counts = result.counts()
    assert counts["H-bonds"] == 2
    assert counts["Hydrophobic interaction"] == 6


def test_every_interaction_column_is_present(result):
    for column in INTERACTION_COLUMNS:
        assert column in result.table.columns


def test_the_ligand_is_not_listed_as_its_own_contact(result):
    for column in INTERACTION_COLUMNS:
        assert not any(r.startswith("UNK") for r in result.residues(column))


def test_fingerprint_matches_the_tables(result):
    fp = fingerprint(EXAMPLE)[0]
    assert fp["hbond_count"] == result.counts()["H-bonds"]
    assert fp["hydrophobic_residues"] == result.residues("Hydrophobic interaction")


def test_a_missing_file_is_reported_clearly():
    with pytest.raises(FileNotFoundError):
        analyse("no_such_complex.pdb")


def test_the_diagram_names_the_ligand_atoms_it_bonds_to():
    """The dashed lines are only meaningful if the detector names a ligand atom."""
    from poseprint.diagram import interaction_records
    found = interaction_records(EXAMPLE)
    assert found, "the example carries a ligand"
    hbonds = [r for r in found[0]["records"] if r["type"] == "hydrogen bond"]
    assert len(hbonds) == 2
    assert {r["resname"] for r in hbonds} == {"ASP"}
    assert all(r["ligand_atom"] for r in hbonds), "each one points at a ligand atom"
    assert all(r["distance"] and 1.5 < r["distance"] < 4.0 for r in hbonds)


def test_the_diagram_renders_and_saves(tmp_path):
    import matplotlib.pyplot as plt

    from poseprint import draw_diagram
    out = tmp_path / "diagram.pdf"
    figure = draw_diagram(EXAMPLE, out=out, title="test")
    plt.close(figure)
    assert out.is_file() and out.stat().st_size > 5000


def test_residue_classes_cover_the_twenty():
    from poseprint.diagram import RESIDUE_CLASS, residue_class
    twenty = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
              "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    assert sorted(RESIDUE_CLASS) == sorted(twenty)
    assert residue_class("HOH") == "other"
