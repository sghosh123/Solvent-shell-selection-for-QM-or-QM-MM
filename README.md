# QM/MM solvent-shell selection

This repository now preserves and ports the snapshot-based workflow in `source/qm_prep_mult_frame_exp_solv.tcl`. The reimplementation is intended to keep the original **selection algorithm unchanged**, not replace it with a generic nearest-water or Cu-centered cutoff rule.

## Original Tcl workflow

For each snapshot, the source script counts unique water residue IDs within 10 Å of the `CU1P DTT DTBD` solute. Its target is 50 waters. If the count is not 50, it retains waters within 3.5 Å of atoms named `S1 S2 S3 S4 N2 O1 O2 CU1P`, then fills the remainder from the other waters in the 10 Å shell, ranked by distance from the geometric center of all solute atoms to the water `OH2` atom. The source's resid-only behavior is preserved for fidelity. Full details and caveats are in `docs/SOURCE_AUDIT.md`.

## Run the Tcl-algorithm port on one snapshot

Install in editable mode:

```bash
python -m pip install -e ".[dev]"
pytest
```

Then run:

```bash
select-solvent-shell tcl-workflow snapshot.pdb \
  --target-waters 50 \
  --output-pdb qm_system_snapshot.pdb \
  --report-json selection_report.json
```

Or without installation:

```bash
PYTHONPATH=src python -m solvent_shell.cli tcl-workflow snapshot.pdb \
  --target-waters 50 --output-pdb qm_system_snapshot.pdb \
  --report-json selection_report.json
```

This command processes a **single PDB snapshot**. Run it separately for each input snapshot; it does not automatically iterate over a trajectory or launch NAMD/ORCA.

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

Regression tests use synthetic PDBs and check the algorithm's main branches. Real-structure parity against VMD remains unverified because the matching original snapshot PDBs were not supplied with the Tcl file. Inspect the JSON report and compare selected residue IDs against VMD before production use.
