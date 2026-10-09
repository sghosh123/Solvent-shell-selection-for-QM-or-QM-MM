# QM/MM solvent-shell selection

This repository preserves and ports the snapshot-based workflow in `source/qm_prep_mult_frame_exp_solv.tcl.

## Original Tcl workflow

For each snapshot, the source script counts unique water residue IDs within 10 Å of the `CU1P DTT DTBD` solute. Its target is 50 waters. If the count is not 50, it retains waters within 3.5 Å of atoms named `S1 S2 S3 S4 N2 O1 O2 CU1P`, then fills the remainder from the other waters in the 10 Å shell, ranked by distance from the geometric center of all solute atoms to the water `OH2` atom. The source's resid-only behavior is preserved for fidelity. Full details and caveats are in `docs/SOURCE_AUDIT.md`.

## Run on one or multiple PDB snapshots

Install the package and tests:

```bash
python -m pip install -e ".[dev]"
pytest
```

One PDB:

```bash
select-solvent-shell tcl-workflow --pdb snapshot.pdb --output-dir qm_frames
```

Multiple PDBs or a wildcard (quote the wildcard so Python expands it):

```bash
select-solvent-shell tcl-workflow --pdb "snapshots/*.pdb" --output-dir qm_frames
```

Each input snapshot produces a separate `<input-stem>_qm_frame.pdb`. A `selection_reports.json` file records the selected residue IDs and counts for every input.

## Run on DCD or XTC trajectories

Install optional trajectory support:

```bash
python -m pip install -e ".[trajectory]"
```

Use a topology that matches the trajectory atom order (PDB, PSF, GRO, etc. supported by MDAnalysis):

```bash
select-solvent-shell tcl-workflow \
  --topology system.pdb \
  --trajectory production.xtc \
  --output-dir qm_frames
```

For DCD, replace `production.xtc` with your `.dcd` file. Each trajectory frame is written as its own PDB, named `frame_000000_qm.pdb`, `frame_000001_qm.pdb`, and so on. Frame numbers are zero-based. To process a subset, use `--start 100 --stop 200 --step 2` (stop is exclusive). The report JSON includes the trajectory name and frame index.

The trajectory topology must contain the solute and solvent atoms, with the residue and atom names expected by the original Tcl script. This tool processes coordinates frame-by-frame; it does not alter the underlying trajectory.

## Generic helper retained

The earlier generic center-atom selector remains available as a separate mode and is not the Tcl algorithm:

```bash
select-solvent-shell generic examples/cu_water_shell_demo.pdb \
  --center-resname CU1 --center-atom CU --water-resnames HOH \
  --oxygen-names OH2 --nearest 2
```

## Important scope distinction

This snapshot-preparation script is not the same as NAMD's live solvent-selection mechanism used during QM/MM dynamics. The available run log describes live solvent selection with `TIP3` and a Cu(I)-relative 4 Å cutoff; those settings should not be conflated with the Tcl script's 10 Å outer shell, 3.5 Å local shell, target of 50 waters, and geometric-center ranking.

## Validation limits

Regression tests use synthetic PDBs and check the algorithm's main branches. Inspect the JSON report and compare selected residue IDs against VMD before production use.
