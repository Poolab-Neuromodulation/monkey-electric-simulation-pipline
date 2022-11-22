import transplant
import os
import argparse
import sys
import shutil
import json

def load_arguments():
    parser = argparse.ArgumentParser(description="Electric field simulation for monkeys")
    parser.add_argument("-o", "--output", help="Directory where the output file is to be saved")
    parser.add_argument("-id", "--ID", help="SUBJECT ID", type=str)
    parser.add_argument("-e", "--electrodes", help="json file for electrodes current dictionary", type=str)
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

    if not args.electrodes:
        print("Should specify subject electrodes. See -h for instructions.")
        sys.exit(1)
    
    OUTPUT_DIR = args.output # out_path
    SUBID = args.ID
    json_path = args.electrodes

    loc = {'Pz': 0, 'Cz': 0, 'Fz': 0, 'Fpz': 0, 'O1': 0, 'C3': 0, 'F1': 0, 'O2': 0, 'C4': 0, 'F2': 0, 'P1': 0, 'P2': 0, 'T1': 0, 'T2': 0,}
    ind = {'Pz': 1, 'Cz': 2, 'Fz': 3, 'Fpz': 4, 'O1': 5, 'C3': 6, 'F1': 7, 'O2': 8, 'C4': 9, 'F2': 10, 'P1': 11, 'P2': 12, 'T1': 13, 'T2': 14,}
    with open(json_path, encoding='utf-8') as f:
        dic = json.load(f)

    ind_e = []
    current = []    
    for key in dic.keys():
        ind_e.append(ind[key])
        loc[key] = dic[key]
        current.append(loc[key])

    ind_e_str = '['
    current_str = '['

    for i in range(len(ind_e)):
        ind_e_str += str(ind_e[i]) + ';'
        current_str += str(current[i]) + ';'
    ind_e_str = ind_e_str[:-1] + ']'
    current_str = current_str[:-1] + ']'


    shutil.copy("./conduct.mat", os.path.join(OUTPUT_DIR, SUBID, 'simu'))

    # coor_A
    print("Coord A")
    os.system('fsleyes '+os.path.join(OUTPUT_DIR, SUBID, 'simu', 'label_tissue.nii'))
            
    A_y = input('Coords A y_axis: ')
    A_z = input('Coords A z_axis: ')

    # boundary            
    boundary_x_low = input('boundary x_low: ')
    boundary_x_high = input('boundary x_high: ')
    boundary_y = input('boundary y plane: ')
    boundary_z = input('boundary z plane: ')


    matlab = transplant.Matlab()

    matlab.simu_monkey_TI(os.path.join(OUTPUT_DIR, SUBID, 'simu'), SUBID, A_y, A_z, boundary_x_low, boundary_x_high, boundary_y, boundary_z, ind_e_str, current_str)
