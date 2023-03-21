import simnibs.mesh_tools.mesh_io as mesh_io
import os
import argparse
import sys

def load_arguments():
    parser = argparse.ArgumentParser(description="Electric field simulation for monkeys")
    parser.add_argument("-o", "--output", help="Directory where the output file is to be saved")
    parser.add_argument("-id", "--ID", help="SUBJECT ID", type=str)
    args = parser.parse_args()
    return args




if __name__ == "__main__":
    args = load_arguments()

    if not args.output:
        print("Should specify output. See -h for instructions.")
        sys.exit(1)
    
    if not args.ID:
        print("Should specify subject ID. See -h for instructions.")
        sys.exit(1)
    
    OUTPUT_DIR = args.output # out_path
    SUBID = args.ID

msh = mesh_io.read_msh(os.path.join( OUTPUT_DIR, SUBID, 'simu', 'msh_tmp', f'{SUBID}.msh'))
msh.reconstruct_surfaces()

mesh_io.write_msh(msh, os.path.join( OUTPUT_DIR, SUBID, 'simu', 'msh_tmp', f'{SUBID}_sim.msh'))

