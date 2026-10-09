"""Command-line interface for PDB, multi-PDB, and trajectory processing."""
import argparse
import csv
import json
from pathlib import Path
import glob

from .selection import select_waters, write_pdb, select_tcl_snapshot, select_tcl_atoms, write_atom_lines, Atom


def _pdb_atom_from_mda(atom, serial):
    """Convert an MDAnalysis atom to the Atom record used by the source-faithful selector."""
    name = str(atom.name)
    resname = str(atom.resname)
    resid = str(int(atom.resid))
    chain = str(getattr(atom, "chainID", "") or "")[:1]
    x, y, z = (float(v) for v in atom.position)
    element = str(getattr(atom, "element", "") or "").strip().upper()
    if not element:
        element = "".join(c for c in name if c.isalpha())[:1].upper()
    # Fixed-column PDB format. Preserve common four-character atom names.
    atom_field = f"{name:>4s}" if len(name) >= 4 else f" {name:<3s}"
    line = (f"ATOM  {serial:5d} {atom_field}{resname:>4s}{chain:1s}{int(atom.resid):4d} "
            f"   {x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {element:>2s}  ")
    return Atom(line, serial, name, resname, chain, resid, "", x, y, z, element)


def _trajectory_atoms(universe):
    return [_pdb_atom_from_mda(atom, i + 1) for i, atom in enumerate(universe.atoms)]


def _run_tcl_pdbs(paths, output_dir, target_waters):
    output_dir.mkdir(parents=True, exist_ok=True)
    reports = []
    for index, path in enumerate(paths):
        atoms, report = select_tcl_snapshot(path, target_waters=target_waters)
        stem = Path(path).stem
        out = output_dir / f"{stem}_qm_frame.pdb"
        write_atom_lines(atoms, out)
        report["output_pdb"] = str(out)
        reports.append(report)
        print(f"[{index + 1}/{len(paths)}] {path}: {report['selected_water_count_unique_resid']} waters -> {out}")
    (output_dir / "selection_reports.json").write_text(json.dumps(reports, indent=2) + "\n", encoding="utf-8")


def _run_tcl_trajectory(topology, trajectory, output_dir, target_waters, start, stop, step):
    try:
        import MDAnalysis as mda
    except ImportError as exc:
        raise SystemExit("Trajectory support requires MDAnalysis. Install with: python -m pip install -e '.[trajectory]'") from exc
    universe = mda.Universe(topology, trajectory)
    n_frames = len(universe.trajectory)
    start = 0 if start is None else start
    stop = n_frames if stop is None else min(stop, n_frames)
    if start < 0 or stop < 0 or step < 1 or start >= stop:
        raise SystemExit(f"Invalid frame range start={start}, stop={stop}, step={step}; trajectory has {n_frames} frames")
    output_dir.mkdir(parents=True, exist_ok=True)
    reports = []
    written = 0
    for frame in range(start, stop, step):
        universe.trajectory[frame]
        atoms = _trajectory_atoms(universe)
        selected, report = select_tcl_atoms(atoms, source_label=f"{trajectory}:frame={frame}", target_waters=target_waters)
        out = output_dir / f"frame_{frame:06d}_qm.pdb"
        write_atom_lines(selected, out)
        report.update({"trajectory": str(trajectory), "topology": str(topology), "frame": frame, "output_pdb": str(out)})
        reports.append(report)
        written += 1
        print(f"[{written}] frame {frame}: {report['selected_water_count_unique_resid']} waters -> {out}")
    (output_dir / "selection_reports.json").write_text(json.dumps(reports, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {written} single-frame PDB files to {output_dir}")


def main():
    p = argparse.ArgumentParser(description="Select QM-region waters using the original Tcl algorithm or a generic selector.")
    sub = p.add_subparsers(dest="command", required=True)

    g = sub.add_parser("generic", help="Generic nearest/cutoff selection (legacy mode)")
    g.add_argument("pdb"); g.add_argument("--center-resname", required=True); g.add_argument("--center-atom", required=True)
    g.add_argument("--center-chain"); g.add_argument("--center-resid")
    g.add_argument("--water-resnames", nargs="+", default=["TIP3", "HOH", "WAT"]); g.add_argument("--oxygen-names", nargs="+", default=["OH2", "O", "OW"])
    gg = g.add_mutually_exclusive_group(required=True); gg.add_argument("--nearest", type=int); gg.add_argument("--cutoff", type=float)
    g.add_argument("--output-pdb", default="selected_waters.pdb"); g.add_argument("--report-csv", default="selected_waters.csv")

    t = sub.add_parser("tcl-workflow", help="Apply the Tcl selection algorithm to one or more PDBs or a DCD/XTC trajectory")
    source = t.add_mutually_exclusive_group(required=True)
    source.add_argument("--pdb", nargs="+", help="One or more PDBs; shell wildcards are supported when quoted")
    source.add_argument("--topology", help="Topology PDB/PSF/GRO for a trajectory")
    t.add_argument("--trajectory", help="DCD/XTC (or another MDAnalysis-supported trajectory); requires --topology")
    t.add_argument("--output-dir", default="qm_frames")
    t.add_argument("--target-waters", type=int, default=50)
    t.add_argument("--start", type=int, default=None, help="First trajectory frame (zero-based, inclusive)")
    t.add_argument("--stop", type=int, default=None, help="Stop frame (zero-based, exclusive)")
    t.add_argument("--step", type=int, default=1, help="Process every Nth frame")

    a = p.parse_args()
    if a.command == "tcl-workflow":
        outdir = Path(a.output_dir)
        if a.pdb:
            if a.trajectory:
                raise SystemExit("Use --pdb for PDB inputs or --topology plus --trajectory for trajectories, not both.")
            paths = []
            for item in a.pdb:
                matches = sorted(glob.glob(item))
                paths.extend(matches if matches else ([item] if Path(item).is_file() else []))
            # Avoid duplicate inputs while preserving order.
            paths = list(dict.fromkeys(paths))
            if not paths:
                raise SystemExit("No input PDB files found.")
            missing = [x for x in paths if not Path(x).is_file()]
            if missing:
                raise SystemExit(f"Input PDB not found: {missing[0]}")
            _run_tcl_pdbs(paths, outdir, a.target_waters)
        else:
            if not a.trajectory:
                raise SystemExit("--topology must be accompanied by --trajectory.")
            _run_tcl_trajectory(a.topology, a.trajectory, outdir, a.target_waters, a.start, a.stop, a.step)
    else:
        center, waters = select_waters(a.pdb, a.center_resname, a.center_atom, a.water_resnames, a.oxygen_names,
                                       a.nearest, a.cutoff, a.center_chain, a.center_resid)
        write_pdb(waters, a.output_pdb)
        with open(a.report_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f); w.writerow(["center_resname", "center_atom", "water_resname", "chain", "resid", "oxygen_atom", "distance_angstrom"])
            for x in waters: w.writerow([center.resname, center.name, x.resname, x.chain, x.resid, x.oxygen, f"{x.distance:.4f}"])
        print(f"Selected {len(waters)} waters; wrote {a.output_pdb} and {a.report_csv}")

if __name__ == "__main__":
    main()
