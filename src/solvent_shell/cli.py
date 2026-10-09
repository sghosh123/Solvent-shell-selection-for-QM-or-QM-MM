import argparse, csv, json
from pathlib import Path
from .selection import select_waters, write_pdb, select_tcl_snapshot, write_atom_lines

def main():
    p=argparse.ArgumentParser(description="Select waters from a PDB using either generic selection or the original Tcl workflow.")
    sub=p.add_subparsers(dest="command", required=True)
    g=sub.add_parser("generic", help="Generic nearest/cutoff selection (legacy mode)")
    g.add_argument("pdb"); g.add_argument("--center-resname",required=True); g.add_argument("--center-atom",required=True)
    g.add_argument("--center-chain"); g.add_argument("--center-resid")
    g.add_argument("--water-resnames",nargs="+",default=["TIP3","HOH","WAT"]); g.add_argument("--oxygen-names",nargs="+",default=["OH2","O","OW"])
    gg=g.add_mutually_exclusive_group(required=True); gg.add_argument("--nearest",type=int); gg.add_argument("--cutoff",type=float)
    g.add_argument("--output-pdb",default="selected_waters.pdb"); g.add_argument("--report-csv",default="selected_waters.csv")
    t=sub.add_parser("tcl-workflow", help="Faithful port of qm_prep_mult_frame_exp_solv.tcl selection logic for one snapshot")
    t.add_argument("pdb"); t.add_argument("--target-waters",type=int,default=50); t.add_argument("--output-pdb",default="qm_system.pdb"); t.add_argument("--report-json",default="selection_report.json")
    a=p.parse_args()
    if a.command=="tcl-workflow":
        atoms,report=select_tcl_snapshot(a.pdb,target_waters=a.target_waters)
        write_atom_lines(atoms,a.output_pdb); Path(a.report_json).write_text(json.dumps(report,indent=2)+"\n")
        print(f"Tcl workflow: outer shell={report['outer_water_count']} unique residues; selected={report['selected_water_count_unique_resid']}; wrote {a.output_pdb}")
        if report['outer_water_count']<a.target_waters: print("WARNING: fewer waters than target within 10 Å; original logic does not add waters beyond that shell.")
        if report['selected_water_count_unique_resid']>a.target_waters: print("WARNING: selected count exceeds target; local-shell rule can make this possible in the source algorithm.")
    else:
        center,waters=select_waters(a.pdb,a.center_resname,a.center_atom,a.water_resnames,a.oxygen_names,a.nearest,a.cutoff,a.center_chain,a.center_resid)
        write_pdb(waters,a.output_pdb)
        with open(a.report_csv,"w",newline="",encoding="utf-8") as f:
            w=csv.writer(f); w.writerow(["center_resname","center_atom","water_resname","chain","resid","oxygen_atom","distance_angstrom"])
            for x in waters: w.writerow([center.resname,center.name,x.resname,x.chain,x.resid,x.oxygen,f"{x.distance:.4f}"])
        print(f"Selected {len(waters)} waters; wrote {a.output_pdb} and {a.report_csv}")
if __name__=="__main__": main()
