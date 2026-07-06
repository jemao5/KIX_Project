import sys
import os
from schrodinger import structure
import glob
from schrodinger.structutils import analyze
from schrodinger.structutils.interactions import hbond
from schrodinger.structutils.interactions.pi import find_pi_pi_interactions
import argparse
import pickle

def find_hbond(chain1_st, chain2_st):
    hbonds = hbond.get_hydrogen_bonds(chain1_st, st2=chain2_st,
                                      max_dist=2.8, 
                                      min_donor_angle=120.0, 
                                      min_acceptor_angle=90.0)
    return hbonds

def pi_pi(chain1_st, chain2_st):
    pi_pi = find_pi_pi_interactions(chain1_st, struct2=chain2_st)
    return pi_pi


def find_all_interactions(out_path, chain_id1, chain_id2):
    st = next(structure.StructureReader(out_path))

    # Formulate the ASL query
    # Note: ASL expects double quotes around the chain name
    asl_expr_chain_id1 = f'chain.name "{chain_id1}"'
    asl_expr_chain_id2 = f'chain.name "{chain_id2}"'

    # Get the list of atom indices (returns a list of integers)
    atom_set_1 = analyze.evaluate_asl(st, asl_expr_chain_id1)
    atom_set_2 = analyze.evaluate_asl(st, asl_expr_chain_id2)
    
    chain1_st = st.extract(atom_set_1)
    chain2_st = st.extract(atom_set_2)
    
    hbonds = find_hbond(chain1_st, chain2_st)
    pi_pi_interactions = pi_pi(chain1_st, chain2_st)
    
    hbonds_info = []
    for hb in hbonds:
        donor_atom = hb[0]
        acceptor_atom = hb[1]
        donor_resnum = donor_atom.resnum
        acceptor_resnum = acceptor_atom.resnum
        donor_resname = donor_atom.pdbres.strip()
        acceptor_resname = acceptor_atom.pdbres.strip()
        donor_chain = donor_atom.chain
        acceptor_chain = acceptor_atom.chain
        if donor_chain == chain_id1:
            hbonds_info.append((f"{donor_resname} {donor_resnum} {donor_chain}", f"{acceptor_resname} {acceptor_resnum} {acceptor_chain}"))
        else:
            hbonds_info.append((f"{acceptor_resname} {acceptor_resnum} {acceptor_chain}", f"{donor_resname} {donor_resnum} {donor_chain}"))
    
    pi_pi_info = []
    for pi in pi_pi_interactions:
        ring1_atom_index = pi.ring1.atoms
        ring2_atom_index = pi.ring2.atoms

        r1_atom_idx = ring1_atom_index[0]
        r2_atom_idx = ring2_atom_index[0]

        r1_atom = pi.struct1.atom[r1_atom_idx]
        r2_atom = pi.struct2.atom[r2_atom_idx]
        r1_resnum = r1_atom.resnum
        r2_resnum = r2_atom.resnum
        r1_resname = r1_atom.pdbres.strip()
        r2_resname = r2_atom.pdbres.strip()
        r1_chain = r1_atom.chain
        r2_chain = r2_atom.chain
        if r1_chain == chain_id1:
            pi_pi_info.append((f"{r1_resname} {r1_resnum} {r1_chain}", f"{r2_resname} {r2_resnum} {r2_chain}"))
        else:
            pi_pi_info.append((f"{r2_resname} {r2_resnum} {r2_chain}", f"{r1_resname} {r1_resnum} {r1_chain}"))
    return hbonds_info, pi_pi_info

def parse_args():
    parser = argparse.ArgumentParser(description="find h-bond and pi stacking interactions for a chunk")
    parser.add_argument("list_file", help="file listing structure paths (one per line)")
    parser.add_argument("pickle_out_path", help="output pickle for this chunk")
    parser.add_argument("--start", type=int, required=True, help="1-based start line (inclusive)")
    parser.add_argument("--end", type=int, required=True, help="1-based end line (inclusive)")
    return parser.parse_args()

def main():
    args = parse_args()

    # read the slice of lines [start, end] from the list file
    with open(args.list_file) as f:
        lines = [ln.strip() for ln in f if ln.strip()]


    chunk_paths = lines[args.start - 1 : args.end]

    all_data = []
    for pdb_path in chunk_paths:
        maegz_path = pdb_path[:-4] + "_prepared.maegz"   # strip ".pdb", add "_prepared.maegz"
        try:
            hb_info, pi_pi_info = find_all_interactions(maegz_path, "A", "B")
            all_data.append({"hbond": hb_info, "pi_pi": pi_pi_info, "file": maegz_path})
        except Exception as e:
            print(f"FAILED {maegz_path}: {e}", file=sys.stderr)
            continue

    with open(args.pickle_out_path, "wb") as f:
        pickle.dump(all_data, f)

if __name__ == "__main__":
    main()
