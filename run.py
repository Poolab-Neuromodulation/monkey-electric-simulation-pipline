import os
import argparse
import sys
import shutil

def load_arguments():
    parser = argparse.ArgumentParser(description="PREEMACS : Pipeline for preprocessing and extraction of the macaque brain surface")
    parser.add_argument("-t1", "--t1_path", help="Path to T1w file(s)")
    parser.add_argument("-ct", "--ct_path", help="Path to CT file(s)")
    parser.add_argument("-o", "--output", help="Directory where the output file is to be saved")
    parser.add_argument("-id", "--ID", help="SUBJECT ID", type=str)
    parser.add_argument('-sphinx', "--Sphinx", help='Reorient to sphinx position', action='store_true')
    parser.add_argument('-mc', "--Manual_Crop", help='if the automatic crop fails (default) do manual crop. Open fslview to id the coordinates (x y z) of anterior and posterior commisures', action='store_true')
    parser.add_argument('-mcc', "--Manual_Crop_coord", help='if the automatic crop fails add coordinates (x y z) of anterior and posterior commisure. After perform all module 1 you can find the file in out_path/file_to_coords.txt. You must add the coords here and run again using this option.', action='store_true')
    parser.add_argument('-fmc', "--Fix_Manual_Crop", help='Manually adjust the cropping range of (x y z) coordinates', action='store_true')
    parser.add_argument('-av_FS', "--av_FS", help='FS method or HCP method (default)', action='store_true')
    parser.add_argument('-asym', "--asym", help='asym template', action='store_true')
    parser.add_argument('-qc_LR', "--qc_LR", help='Based on vitame E capsule verify left-rigth side', action='store_true')
    parser.add_argument('-tmp', "--TMP", help='do not remove temporal files', action='store_true')
    parser.add_argument("-e", "--electrodes", help="json file for electrodes current dictionary", type=str)
    parser.add_argument('-ti', "--TI", help='Reorient to sphinx position', action='store_true')
    args = parser.parse_args()
    return args
    
if __name__ == "__main__":
    args = load_arguments()

    if not args.t1_path:
        print("Should specify t1_path. See -h for instructions.")
        sys.exit(1)

    if not args.ct_path:
        print("Should specify ct_path. See -h for instructions.")
        sys.exit(1)

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
    T1_image_path = args.t1_path
    CT_image_path = args.ct_path
    json_path = args.electrodes    
        
    if not args.Fix_Manual_Crop:   
        os.system(f'python preprocess.py -t1 {T1_image_path} -ct {CT_image_path} -id {SUBID} -o {OUTPUT_DIR}') 
    if args.Fix_Manual_Crop:
        os.system(f'python preprocess.py -t1 {T1_image_path} -ct {CT_image_path} -id {SUBID} -o {OUTPUT_DIR} -fmc')    
    os.system(f'python ct_segmentation.py -id {SUBID} -o {OUTPUT_DIR}')

    if not args.TI:
        os.system(f'python simu_monkey.py -id {SUBID} -o {OUTPUT_DIR} -e {json_path}')

    if args.TI:
        os.system(f'python simu_monkey_TI.py -id {SUBID} -o {OUTPUT_DIR}')
