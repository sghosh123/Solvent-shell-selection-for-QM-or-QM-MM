from pathlib import Path
import pytest
from solvent_shell.selection import select_waters, write_pdb

FIXTURE = Path(__file__).parent / "fixtures" / "cu_water_shell_demo.pdb"

def test_nearest_two_waters():
    center, waters = select_waters(FIXTURE, "CU1", "CU", ["HOH"], ["OH2"], nearest=2)
    assert center.resid == "900"
    assert [w.resid for w in waters] == ["101", "102"]
    assert [round(w.distance, 1) for w in waters] == [2.0, 3.2]
    assert all(len(w.atoms) == 3 for w in waters)

def test_cutoff_four_angstrom():
    _, waters = select_waters(FIXTURE, "CU1", "CU", ["HOH"], ["OH2"], cutoff=4.0)
    assert [w.resid for w in waters] == ["101", "102"]

def test_modes_mutually_exclusive():
    with pytest.raises(ValueError):
        select_waters(FIXTURE, "CU1", "CU", nearest=2, cutoff=4.0)

def test_complete_residue_output(tmp_path):
    _, waters = select_waters(FIXTURE, "CU1", "CU", ["HOH"], ["OH2"], nearest=1)
    out = tmp_path / "selected.pdb"
    write_pdb(waters, out)
    s = out.read_text()
    assert "OH2  HOH" in s and "H1   HOH" in s and "H2   HOH" in s
    assert s.rstrip().endswith("END")

from solvent_shell.selection import select_tcl_snapshot

def _pdb_atom(serial, name, resname, resid, x, y, z, element, chain='A'):
    return (f"ATOM  {serial:5d} {name:>4s} {resname:>4s}{chain}{resid:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {element:>2s}  ")

def _make_tcl_fixture(path, nwaters, local_count=1):
    lines=[]; serial=1
    # Solute geometry centered near origin; one local atom is the shell anchor.
    for name,xyz in [('CU1P',(0,0,0)),('S1',(0,0,0)),('D1',(1,0,0))]:
        lines.append(_pdb_atom(serial,name,'CU1P',900,*xyz,'C' if name=='D1' else 'CU')); serial+=1
    # Place one water close to S1; remaining waters within outer 10 Å, ordered by center distance.
    for i in range(nwaters):
        x=1.0+i*0.1
        resid=101+i
        for name,dx,element in [('OH2',0.0,'O'),('H1',0.1,'H'),('H2',-0.1,'H')]:
            lines.append(_pdb_atom(serial,name,'TIP3',resid,x+dx,0,0,element)); serial+=1
    path.write_text('\n'.join(lines)+'\nEND\n')

def test_tcl_algorithm_exact_target_keeps_outer_shell(tmp_path):
    pdb=tmp_path/'snapshot.pdb'; _make_tcl_fixture(pdb,50)
    atoms, report=select_tcl_snapshot(pdb)
    assert report['outer_water_count']==50
    assert report['selected_water_count_unique_resid']==50
    assert report['mode']=='exact_target_keep_outer_shell'
    assert sum(a.resname=='TIP3' and a.name=='OH2' for a in atoms)==50

def test_tcl_algorithm_local_shell_then_nearest_fill(tmp_path):
    pdb=tmp_path/'snapshot.pdb'; _make_tcl_fixture(pdb,55)
    atoms, report=select_tcl_snapshot(pdb)
    assert report['outer_water_count']==55
    assert report['selected_water_count_unique_resid']==50
    assert report['mode']=='local_shell_plus_center_ranked_outer_waters'
    assert sum(a.resname=='TIP3' and a.name=='OH2' for a in atoms)==50
    # The local-shell water is retained even if not among the center-nearest candidates.
    assert '101' in report['selected_water_resids']

def test_tcl_algorithm_does_not_invent_waters_when_outer_count_below_target(tmp_path):
    pdb=tmp_path/'snapshot.pdb'; _make_tcl_fixture(pdb,8)
    atoms, report=select_tcl_snapshot(pdb)
    assert report['outer_water_count']==8
    assert report['selected_water_count_unique_resid']==8
    assert sum(a.resname=='TIP3' and a.name=='OH2' for a in atoms)==8


def test_tcl_atom_api_matches_pdb_api():
    from solvent_shell.selection import read_atoms, select_tcl_atoms
    atoms_file, report_file = select_tcl_snapshot(FIXTURE, target_waters=2, solute_resnames=("CU1",), water_resnames=("HOH",))
    atoms_memory, report_memory = select_tcl_atoms(read_atoms(FIXTURE), source_label="fixture", target_waters=2, solute_resnames=("CU1",), water_resnames=("HOH",))
    assert [(a.name, a.resname, a.resid, a.x, a.y, a.z) for a in atoms_file] == [
        (a.name, a.resname, a.resid, a.x, a.y, a.z) for a in atoms_memory
    ]
    assert report_file["selected_water_resids"] == report_memory["selected_water_resids"]


def test_mult_pdb_cli_writes_one_output_per_input(tmp_path):
    from solvent_shell.cli import _run_tcl_pdbs
    first = tmp_path / "one.pdb"
    second = tmp_path / "two.pdb"
    _make_tcl_fixture(first, 50)
    _make_tcl_fixture(second, 55)
    outdir = tmp_path / "out"
    _run_tcl_pdbs([str(first), str(second)], outdir, 50)
    assert (outdir / "one_qm_frame.pdb").is_file()
    assert (outdir / "two_qm_frame.pdb").is_file()
    reports = __import__("json").loads((outdir / "selection_reports.json").read_text())
    assert len(reports) == 2
    assert [r["outer_water_count"] for r in reports] == [50, 55]


def test_trajectory_adapter_preserves_atom_fields():
    pytest.importorskip("MDAnalysis")
    import MDAnalysis as mda
    from solvent_shell.cli import _pdb_atom_from_mda
    u = mda.Universe(str(FIXTURE))
    atom = _pdb_atom_from_mda(u.atoms[0], 1)
    assert atom.name == u.atoms[0].name
    assert atom.resname == u.atoms[0].resname
    assert atom.resid == str(u.atoms[0].resid)
    assert abs(atom.x - float(u.atoms[0].position[0])) < 1e-3
