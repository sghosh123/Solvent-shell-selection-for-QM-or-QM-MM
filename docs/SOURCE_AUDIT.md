# Source audit and algorithm provenance

## Original source preserved
`source/qm_prep_mult_frame_exp_solv.tcl` is the user's original Tcl/VMD script, copied without modification. It is the primary source for the batch snapshot-selection workflow implemented here.

## Tcl algorithm, as implemented
For each snapshot, the original script:

1. Counts unique water `resid` values for waters within 10 Å of any atom in `resname CU1P DTT DTBD`.
2. Uses a fixed target of 50 waters.
3. If the count is exactly 50, writes the solute plus all waters in that 10 Å selection.
4. Otherwise, selects a local shell of waters within 3.5 Å of solute atoms named `S1 S2 S3 S4 N2 O1 O2 CU1P`.
5. Removes local-shell residue IDs from the outer-shell candidates.
6. Calculates the geometric center of all solute atoms and ranks remaining candidate waters by the distance from that center to each water's `OH2` atom.
7. Retains the nearest candidates needed to reduce the outer-shell count by `count - 50`, then combines them with the local shell and writes the solute plus selected waters.

The Python port intentionally preserves these selection rules, the 10 Å and 3.5 Å cutoffs, the target of 50, the geometric-center/`OH2` ranking, and the original resid-only selection behavior. It does not substitute a Cu-centered cutoff or a generic nearest-water algorithm.

## Source issues preserved or made explicit

- The Tcl script loops over snapshots 0–103, then separately loads snapshot 104 for counting. This is an apparent indexing inconsistency: the count list includes 105 entries, but later processing only iterates over 0–103. The Python CLI processes one explicitly supplied snapshot; a future batch wrapper should default to the intended 0–103 range and flag the extra count pass rather than silently reproducing the mismatch.
- If fewer than 50 waters are within 10 Å, the Tcl list operation retains all candidate waters and does not add waters from outside 10 Å. The port preserves that behavior and reports the count rather than inventing extra waters.
- If the local shell itself exceeds 50 unique residue IDs, the algorithm may return more than 50. The port preserves the local-shell-first rule and reports the actual count.
- The original uses residue ID alone for de-duplication and final selection, without chain qualification. The port preserves this behavior for fidelity; structures with repeated residue IDs across chains should be checked carefully.
- The original creates selections without consistently deleting them, and the expression around `$solv_coord` at the `veclength` call is fragile Tcl syntax. These are implementation hazards, not algorithmic selection rules; the Python implementation avoids them.

## Validation scope
Synthetic regression tests cover the exact-target branch, the local-shell-plus-nearest-fill branch, and the below-target behavior. These validate the Python implementation against the logic expressed in the Tcl source, not against VMD output from real snapshots. No original PDB snapshots were available alongside the Tcl file, so a structure-by-structure parity check against VMD remains outstanding.

## Separate QM/MM live-solvent mechanism
Other available run-log/manuscript material describes NAMD/ORCA live solvent selection, including a Cu(I)-relative 4 Å cutoff. That is a distinct workflow. This repository's Tcl snapshot preprocessing algorithm is not represented as equivalent to the NAMD live-solvent selection feature.
