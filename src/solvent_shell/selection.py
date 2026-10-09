"""PDB parsing and water selection; includes a faithful Tcl-workflow port."""
from dataclasses import dataclass
from pathlib import Path
import math

@dataclass(frozen=True)
class Atom:
    line: str; serial: int; name: str; resname: str; chain: str; resid: str; icode: str
    x: float; y: float; z: float; element: str
    @property
    def key(self): return (self.chain, self.resid, self.icode, self.resname)

@dataclass(frozen=True)
class Water:
    chain: str; resid: str; resname: str; oxygen: str; distance: float; atoms: tuple

def parse_atom(line):
    if not line.startswith(("ATOM  ", "HETATM")): return None
    try:
        name=line[12:16].strip(); element=line[76:78].strip().upper() if len(line)>=78 else ""
        if not element: element="".join(c for c in name if c.isalpha())[:1].upper()
        return Atom(line.rstrip("\n"), int(line[6:11]), name, line[17:21].strip(), line[21:22].strip(),
                    line[22:26].strip(), line[26:27].strip(), float(line[30:38]), float(line[38:46]),
                    float(line[46:54]), element)
    except (ValueError, IndexError): return None

def read_atoms(path):
    atoms=[]; model_seen=False
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("MODEL"):
                if model_seen: break
                model_seen=True; continue
            if line.startswith("ENDMDL") and model_seen: break
            a=parse_atom(line)
            if a: atoms.append(a)
    return atoms

def group_residues(atoms):
    groups={}
    for a in atoms: groups.setdefault(a.key, []).append(a)
    return groups

def distance(a,b): return math.sqrt((a.x-b.x)**2+(a.y-b.y)**2+(a.z-b.z)**2)
def geom_center(atoms):
    if not atoms: raise ValueError("Cannot calculate geometric center of an empty selection")
    n=len(atoms); return (sum(a.x for a in atoms)/n, sum(a.y for a in atoms)/n, sum(a.z for a in atoms)/n)
def point_distance(point, atom): return math.sqrt(sum((x-y)**2 for x,y in zip(point,(atom.x,atom.y,atom.z))))
def residue_id_set(groups):
    # The source Tcl uses resid alone, not chain/resname-qualified residue identity.
    return {a.resid for group in groups.values() for a in group}
def write_atom_lines(atoms, path):
    lines=[]
    for serial,a in enumerate(atoms,1): lines.append(f"{a.line[:6]}{serial:5d}{a.line[11:]}")
    Path(path).write_text("\n".join(lines)+("\n" if lines else "")+"END\n", encoding="utf-8")

def select_tcl_snapshot(pdb_path, target_waters=50, solute_resnames=("CU1P","DTT","DTBD"),
                        local_atom_names=("S1","S2","S3","S4","N2","O1","O2","CU1P"),
                        water_resnames=("TIP3","HOH","WAT"), oxygen_name="OH2",
                        outer_cutoff=10.0, local_cutoff=3.5):
    """Port the Tcl selection logic without changing its selection algorithm.

    Counts waters within outer_cutoff of any solute atom. If that count equals
    target_waters, retains all such waters. Otherwise it retains all waters
    within local_cutoff of named solute atoms, then fills to the target by
    keeping the nearest-to-solute-geometric-center OH2 waters from the outer
    shell candidates. Distances are in Angstrom, as in the Tcl/VMD script.
    The original uses resid alone when writing the final selection; this port
    preserves that behavior.
    """
    atoms=read_atoms(pdb_path); groups=group_residues(atoms)
    solute=[a for a in atoms if a.resname in solute_resnames]
    if not solute: raise ValueError(f"No solute atoms found for resnames {solute_resnames}")
    water_names=set(water_resnames)
    waters={key:group for key,group in groups.items() if group[0].resname in water_names}
    solute_xyz=[(a.x,a.y,a.z) for a in solute]
    local_centers=[a for a in solute if a.name in local_atom_names]
    # VMD 'water and within cutoff of (...)' is atom-based: any atom of a
    # water residue may satisfy the cutoff relative to any selected solute atom.
    def near_group(wgroup, centers, cutoff):
        return any(distance(w,c) <= cutoff for w in wgroup for c in centers)
    outer={key:g for key,g in waters.items() if near_group(g,solute,outer_cutoff)}
    outer_resids={g[0].resid for g in outer.values()}
    count=len(outer_resids)
    if count == target_waters:
        selected_resids=set(outer_resids); mode="exact_target_keep_outer_shell"
    else:
        local={key:g for key,g in waters.items() if local_centers and near_group(g,local_centers,local_cutoff)}
        local_resids={g[0].resid for g in local.values()}
        candidates={key:g for key,g in outer.items() if g[0].resid not in local_resids}
        center=geom_center(solute)
        distances={}
        for key,g in candidates.items():
            oxygen=[a for a in g if a.name == oxygen_name]
            # The Tcl script assumes one OH2 atom and would fail/behave oddly otherwise.
            if oxygen: distances[g[0].resid]=point_distance(center,oxygen[0])
        # Tcl sorts ascending and keeps lrange 0 end-(count-target). If count <
        # target, Tcl's oversized end index effectively retains every candidate;
        # preserve this behavior instead of silently inventing farther waters.
        delta=count-target_waters
        sorted_resids=sorted(distances, key=lambda r:(distances[r],r))
        keep_n=max(0, len(sorted_resids)-delta) if delta >= 0 else len(sorted_resids)
        retained=set(sorted_resids[:keep_n])
        selected_resids=local_resids | retained
        mode="local_shell_plus_center_ranked_outer_waters"
    # Match Tcl's resid-only selection, which can include same-resid atoms from
    # another chain/residue. This is intentional fidelity, not recommended redesign.
    output=[a for a in atoms if a.resname in solute_resnames or
            (a.resid in selected_resids)]
    report={"input":str(pdb_path),"mode":mode,"outer_water_count":count,
            "target_water_count":target_waters,"selected_water_resids":sorted(selected_resids, key=lambda x:(int(x) if x.lstrip('-').isdigit() else 10**12,x)),
            "selected_water_count_unique_resid":len(selected_resids),"solute_atom_count":len(solute),
            "local_atom_names":list(local_atom_names),"outer_cutoff_angstrom":outer_cutoff,
            "local_cutoff_angstrom":local_cutoff}
    return output, report

def select_waters(pdb_path, center_resname, center_atom, water_resnames=("TIP3","HOH","WAT"),
                  oxygen_names=("OH2","O","OW"), nearest=None, cutoff=None, center_chain=None, center_resid=None):
    """Generic static selector retained for backward compatibility."""
    if (nearest is None)==(cutoff is None): raise ValueError("Specify exactly one of nearest or cutoff.")
    if nearest is not None and nearest<1: raise ValueError("nearest must be >= 1.")
    if cutoff is not None and cutoff<=0: raise ValueError("cutoff must be > 0 Angstrom.")
    atoms=read_atoms(pdb_path); groups=group_residues(atoms)
    centers=[a for a in atoms if a.resname.upper()==center_resname.upper() and a.name.upper()==center_atom.upper()
             and (center_chain is None or a.chain==center_chain) and (center_resid is None or a.resid==str(center_resid))]
    if len(centers)!=1: raise ValueError(f"Expected one center atom; found {len(centers)}. Add center chain/resid if needed.")
    c=centers[0]; found=[]
    for group in groups.values():
        if group[0].resname.upper() not in {x.upper() for x in water_resnames}: continue
        oxy=[a for a in group if a.name.upper() in {x.upper() for x in oxygen_names}]
        if len(oxy)!=1: continue
        d=distance(c,oxy[0]); found.append(Water(oxy[0].chain,oxy[0].resid,oxy[0].resname,oxy[0].name,d,tuple(group)))
    found.sort(key=lambda x:(x.distance,x.chain,x.resid)); chosen=found[:nearest] if nearest is not None else [x for x in found if x.distance<=cutoff]
    return c, chosen

def write_pdb(waters, path):
    # Backwards-compatible writer for generic selector tuples.
    atoms=[]
    for item in waters: atoms.extend(item.atoms)
    write_atom_lines(atoms,path)
