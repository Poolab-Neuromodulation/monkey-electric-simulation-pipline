from curses.ascii import SUB
from importlib.resources import path
import os, time, sys, shutil, re
from telnetlib import TM
from nipype.interfaces import fsl
import nipype.interfaces.mrtrix as mrt
import nipype.interfaces.mrtrix3 as mrt3
from nipype.interfaces.freesurfer import MRIConvert
from nipype.interfaces.ants import N4BiasFieldCorrection
from nipype.interfaces.fsl import ExtractROI
from nipype.testing import anatfile
import copy
import transplant
import nibabel as nib
import numpy as np
from deepbrain import Extractor
import argparse



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
    args = parser.parse_args()
    return args


def run_Brain_Mask(img_path, output_dir, p=0.5):
    img = nib.load(img_path)

    affine = img.affine
    img = img.get_fdata()

    extractor = Extractor()

    now = time.time()
    prob = extractor.run(img)
    print("Extraction time: {0:.2f} secs.".format(time.time() - now))
    mask = prob > p
    brain_mask = (1 * mask).astype(np.uint8)
    brain_mask = nib.Nifti1Image(brain_mask, affine)
    nib.save(brain_mask, os.path.join(output_dir, "brain_mask.nii"))

    brain = img[:]
    brain[~mask] = 0
    brain = nib.Nifti1Image(brain, affine)
    nib.save(brain, os.path.join(output_dir, "brain.nii"))


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

    # PATH
    OUTPUT_DIR = args.output # out_path
    SUBID = args.ID
    T1_image_path = args.t1_path
    CT_image_path = args.ct_path
    templates_path = os.path.join(os.getcwd(), 'templates')

    print("output DIR:", OUTPUT_DIR)
    print('SUBID:', SUBID)
    print('T1_DIR:', T1_image_path)
    print('CT_DIR:', CT_image_path)

    # Create output files
    os.makedirs(os.path.join(OUTPUT_DIR, SUBID), mode=0o777, exist_ok=True)
    os.makedirs(os.path.join(OUTPUT_DIR, SUBID, 'crop'), mode=0o777, exist_ok=True)
    os.makedirs(os.path.join(OUTPUT_DIR, SUBID, 'crop', 'antsREg'), mode=0o777, exist_ok=True)
    os.makedirs(os.path.join(OUTPUT_DIR, SUBID, 'image_conform'), mode=0o777, exist_ok=True)
    os.makedirs(os.path.join(OUTPUT_DIR, SUBID, 'N4_T1'), mode=0o777, exist_ok=True)
    os.makedirs(os.path.join(OUTPUT_DIR, SUBID, 'orig'), mode=0o777, exist_ok=True)
    os.makedirs(os.path.join(OUTPUT_DIR, SUBID, 'reorient'), mode=0o777, exist_ok=True)
    os.makedirs(os.path.join(OUTPUT_DIR, SUBID, 'mask'), mode=0o777, exist_ok=True)
    os.makedirs(os.path.join(OUTPUT_DIR, SUBID, 'tmp'), mode=0o777, exist_ok=True)
    os.makedirs(os.path.join(OUTPUT_DIR, SUBID, 'simu'), mode=0o777, exist_ok=True)

    path_job = os.path.join(OUTPUT_DIR, SUBID)
    reorient_path = os.path.join(OUTPUT_DIR, SUBID, 'reorient')
    TMP = os.path.join(path_job, 'tmp')
    image_conform = os.path.join(OUTPUT_DIR, SUBID, 'image_conform')
    path_crop = os.path.join(OUTPUT_DIR, SUBID, 'crop')
    path_ants_reg = os.path.join(OUTPUT_DIR, SUBID, 'crop', 'antsREg')
    N4_T1_path = os.path.join(OUTPUT_DIR, SUBID, 'N4_T1')
    FSL_DIR = os.environ.get('FSLDIR')

    # file_list = [x for x in os.listdir(T1_image_path) if 'nii.gz' in x]
    # num_ima = 1
    # for filename in file_list:
    try:
        shutil.copyfile(T1_image_path,
                        os.path.join(os.path.join(OUTPUT_DIR, SUBID, 'orig'), 'raw_T1.nii.gz')
                        )
            
    except Exception as e:
        print('error')

    # file_list = [x for x in os.listdir(CT_image_path) if 'nii.gz' in x]
    # for filename in file_list:
    try:
        shutil.copyfile(CT_image_path,
                        os.path.join(os.path.join(OUTPUT_DIR, SUBID, 'simu'), 'CT.nii.gz')
                        )
    except Exception as e:
        print('error')
        


    t0 = time.time()
    print('\n')
    print("====================== runnig Module 1 ======================")
    print('\n')

    # optional: sphinx position
    if args.Sphinx:
        file_list = [x for x in os.listdir(os.path.join(OUTPUT_DIR, SUBID, 'orig')) if 'nii.gz' in x]
        for filename in file_list:
            print(os.path.join(OUTPUT_DIR, SUBID, 'orig',filename))
            cvt = MRIConvert()
            cvt.inputs.in_file = os.path.join(OUTPUT_DIR, SUBID, 'orig',filename)
            cvt.inputs.out_file = os.path.join(OUTPUT_DIR, SUBID, 'orig',filename)
            cvt.inputs.sphinx = True
            cvt.run()


    print("======================================================")
    print("=                                                    =")
    print("=                Volume Orientation                  =")
    print("=                                                    =")
    print("======================================================")
    file_list = [x for x in os.listdir(os.path.join(OUTPUT_DIR, SUBID, 'orig')) if 'nii.gz' in x]
    for filename in file_list:
        reorient = fsl.Reorient2Std(out_file=os.path.join(reorient_path, filename[0:-7]+'_REO.nii.gz'))
        reorient.inputs.in_file = os.path.join(OUTPUT_DIR, SUBID, 'orig', filename)
        res = reorient.run()

        os.system('mrinfo '+os.path.join(reorient_path, filename[0:-7]+'_REO.nii.gz')+' > '+os.path.join(TMP, filename[0:-7]+'_REO.txt'))
        os.system('grep strides '+os.path.join(TMP, filename[0:-7]+'_REO.txt')+' > '+os.path.join(TMP, filename[0:-7]+'_strides.txt'))
        orient = os.popen("awk \'{ print $4 }\' "+os.path.join(TMP, filename[0:-7]+'_strides.txt'))
        orient = int(orient.read())
        correct_orient = 1
        
        if orient != correct_orient:
            print("Wrong orientation...........")
            # fix orientation

            os.system(f'fslorient -deleteorient {reorient_path}/{filename[0:-7]}_REO.nii.gz && fslswapdim {reorient_path}/{filename[0:-7]}_REO.nii.gz  -x y z {TMP}/{filename[0:-7]}_fix_REO.nii.gz')
            os.remove(os.path.join(reorient_path, filename[0:-7]+'_REO.nii.gz'))
            os.system(f'fslorient -setqformcode 1  {TMP}/{filename[0:-7]}_fix_REO.nii.gz')
            try:
                shutil.copyfile(os.path.join(TMP, filename[0:-7]+'_fix_REO.nii.gz'),
                                os.path.join(reorient_path, filename[0:-7]+'_REO.nii.gz')
                                )
            except Exception as e:
                print('error')

        if args.qc_LR:
            if orient != correct_orient:
                print("Wrong orientation...........")
                os.system(f'fslview {reorient_path}/{filename[0:-7]}_REO.nii.gz')
                try:
                    shutil.copyfile(os.path.join(reorient_path, filename[0:-7]+'_REO.nii.gz'),
                                os.path.join(TMP, filename[0:-7]+'_REO_original.nii.gz')
                                )
                except Exception as e:
                    print('error')
                print("> Is the vitamin E capsule on the right side? Y or N or NI (no info) and close fsl")
                orientation = input("Orientation: ")
                if orientation == 'N':
                    os.system(f'fslorient -deleteorient {reorient_path}/{filename[0:-7]}_REO.nii.gz && fslswapdim {reorient_path}/{filename[0:-7]}_REO.nii.gz  x y z {TMP}/{filename[0:-7]}_fix_REO.nii.gz')
                    os.remove(os.path.join(reorient_path, filename[0:-7]+'_REO.nii.gz'))
                    os.system(f'fslorient -setqformcode 1  {TMP}/{filename[0:-7]}_fix_REO.nii.gz')
                    try:
                        shutil.copyfile(os.path.join(TMP, filename[0:-7]+'_fix_REO.nii.gz'),
                                os.path.join(reorient_path, filename[0:-7]+'_REO.nii.gz')
                                )
                    except Exception as e:
                        print('error')
                elif orientation == 'Y':
                    os.system(f'fslorient -deleteorient {reorient_path}/{filename[0:-7]}_REO.nii.gz && fslswapdim {reorient_path}/{filename[0:-7]}_REO.nii.gz  -x y z {TMP}/{filename[0:-7]}_fix_REO.nii.gz')
                    os.remove(os.path.join(reorient_path, filename[0:-7]+'_REO.nii.gz'))
                    os.system(f'fslorient -setqformcode 1  {TMP}/{filename[0:-7]}_fix_REO.nii.gz')
                    try:
                        shutil.copyfile(os.path.join(TMP, filename[0:-7]+'_fix_REO.nii.gz'),
                                os.path.join(reorient_path, filename[0:-7]+'_REO.nii.gz')
                                )
                    except Exception as e:
                        print('error')
                elif orientation == 'NI':
                    os.system(f'fslorient -deleteorient {reorient_path}/{filename[0:-7]}_REO.nii.gz && fslswapdim {reorient_path}/{filename[0:-7]}_REO.nii.gz  x y z {TMP}/{filename[0:-7]}_fix_REO.nii.gz')
                    os.remove(os.path.join(reorient_path, filename[0:-7]+'_REO.nii.gz'))
                    os.system(f'fslorient -setqformcode 1  {TMP}/{filename[0:-7]}_fix_REO.nii.gz')
                    try:
                        shutil.copyfile(os.path.join(TMP, filename[0:-7]+'_fix_REO.nii.gz'),
                                os.path.join(reorient_path, filename[0:-7]+'_REO.nii.gz')
                                )
                    except Exception as e:
                        print('error')


    print('\n')
    print("done...")
    print('\n')
    print('\n')
    print("======================================================")
    print("=                                                    =")
    print("=                  Image Crop                        =")
    print("=                                                    =")
    print("======================================================")
    print('\n')
    # conform
    file = open(os.path.join(image_conform, 'file_to_coords.txt'), 'w')
    matlab = transplant.Matlab()
    file_list = [x for x in os.listdir(reorient_path) if 'nii.gz' in x]
    for filename in file_list:
        nii_in = filename
        nii_out = filename[0:-11]+'_CROP_REO.nii.gz'
        path_out = image_conform
        conform = matlab.conform2(reorient_path+'/', nii_in, image_conform+'/', nii_out)
        file.write(nii_out+'\n')
    file.close()

    T1_reo = nib.load(os.path.join(image_conform, 'raw_T1_CROP_REO.nii.gz'))
    T1_reo_data = T1_reo.get_fdata()




    if not args.Manual_Crop and not args.Manual_Crop_coord and not args.Fix_Manual_Crop:
    # auto crop
        file_list = [x for x in os.listdir(image_conform) if 'T1_CROP_REO.nii.gz' in x]
        for filename in file_list:
            os.system(f'antsRegistrationSyN.sh -d 3 -f {templates_path}/NMT_05.nii.gz -t r -m {image_conform}/{nii_out} -o {path_ants_reg}/{filename[0:-7]}_REG_')

        file_list = [x for x in os.listdir(path_ants_reg) if '_CROP_REO_REG_0GenericAffine.mat' in x]  
        for filename in file_list:
            os.remove(os.path.join(path_ants_reg, filename))
        file_list = [x for x in os.listdir(path_ants_reg) if '_CROP_REO_REG_Warped.nii.gz' in x]   
        for filename in file_list:
            os.remove(os.path.join(path_ants_reg, filename))  

        file_list = [x for x in os.listdir(path_ants_reg) if 'nii.gz' in x]   
        for filename in file_list:
            reg_file = filename
            original_file= os.path.join(image_conform, reg_file[0:-25]+'.nii.gz')
            nii_out1 = os.path.join(TMP, reg_file[0:-7]+'_ones.nii.gz')
            nii_out2 = reg_file[0:-7]+'_square_ones.nii.gz'
            nii_out3 = reg_file[0:-7]+'_crop.nii.gz'
            out_not_square_FOV = reg_file[0:-7]+'_not_square_FOV.nii.gz'
            out_square_FOV = reg_file[0:-7]+'_square_FOV.nii.gz'

            ones_reg_ants = matlab.crop_ants(path_ants_reg+'/', reg_file, nii_out1)
            
            maths = fsl.MultiImageMaths()
            maths.inputs.in_file = nii_out1
            maths.inputs.op_string = "-mul %s"
            maths.inputs.operand_files = original_file
            maths.inputs.out_file = os.path.join(TMP, out_not_square_FOV)
            maths.cmdline
            maths.run()

            newMat = matlab.preemacs_square_crop(TMP+'/',out_not_square_FOV, nii_out2)
            
            maths = fsl.MultiImageMaths()
            maths.inputs.in_file = os.path.join(TMP, nii_out2)
            maths.inputs.op_string = "-mul %s"
            maths.inputs.operand_files = original_file
            maths.inputs.out_file = os.path.join(TMP, out_square_FOV)
            maths.run()

            newMat = matlab.preemacs_autocrop(TMP+'/',out_square_FOV, nii_out3)

            try:
                shutil.copyfile(os.path.join(TMP, nii_out3),
                                os.path.join(path_crop, nii_out3)
                                )
            except Exception as e:
                print('error')

            size = filename[0:-7]+'_crop.txt'
            nii_prefinal_crop = filename[0:-7]+'_crop.nii.gz'
            image_crop = matlab.crop_only_brain_3(TMP+'/', nii_prefinal_crop, size)
            try:
                f = open(os.path.join(TMP, size), 'r')
                var = f.read()
                f.close()
                print(var)
                coord = re.findall(r'\d+', var)



                final_crop = filename[0:-34]+'_final_crop.nii.gz'

                os.system('mrgrid '+var.strip('\n')+' '+os.path.join(TMP, nii_prefinal_crop)+' crop '+os.path.join(TMP, final_crop))

                try:
                    shutil.copyfile(os.path.join(TMP, final_crop),
                                    os.path.join(path_crop, final_crop)
                                    )
                except Exception as e:
                    print('error')
            except Exception as e:
                print('Auto Crop failed. use -mc or -fmc to execute Manual Crop.')
                sys.exit(1)
    
    if args.Manual_Crop and not args.Manual_Crop_coord and not args.Fix_Manual_Crop:  
        print('\n')
        print("======================================================")
        print("=                                                    =")
        print("=                    Manual Crop                     =")
        print("=                                                    =")
        print("======================================================")
        print('\n')

        file_list = [x for x in os.listdir(image_conform) if 'nii.gz' in x]
        for filename in file_list:
            try:
                shutil.copyfile(os.path.join(image_conform, filename),
                                os.path.join(TMP, filename)
                                )
            except Exception as e:
                print('error')
            
            print("> Coords of ANTERIOR COMMISURE (AC) POSTERIOR COMMISURE (PC) 71 88 147 71 71 146 (see the example .png)")
            os.system('fsleyes '+os.path.join(image_conform, filename))
            
            coords = input('Coords: ')
            size = filename[0:-7]+'_crop.txt'
            image_crop = matlab.manual_crop_brain(TMP+'/', filename, coords, size)
            f = open(os.path.join(TMP, size), 'r')
            var = f.read()
            f.close()
            print(var)
            coord = re.findall(r'\d+', var)
           
            final_crop = filename[0:-16]+'_final_crop.nii.gz'



            os.system('mrgrid '+var.strip('\n')+' '+os.path.join(TMP, filename)+' crop '+os.path.join(TMP, final_crop))

            try:
                shutil.copyfile(os.path.join(TMP, final_crop),
                                os.path.join(path_crop, final_crop)
                                )
            except Exception as e:
                print('error')


    if args.Fix_Manual_Crop and not args.Manual_Crop_coord and not args.Manual_Crop:  
        print('\n')
        print("======================================================")
        print("=                                                    =")
        print("=                    Manual Crop                     =")
        print("=                                                    =")
        print("======================================================")
        print('\n')

        file_list = [x for x in os.listdir(image_conform) if 'nii.gz' in x]
        for filename in file_list:
            try:
                shutil.copyfile(os.path.join(image_conform, filename),
                                os.path.join(TMP, filename)
                                )
            except Exception as e:
                print('error')
            
            print("> Enter the crop range for manual adjustment of (x y z) range.")
            os.system('fsleyes '+os.path.join(image_conform, filename))
            
            shift_i_3 = input('x_start: ')
            x_end = input('x_end: ')
            shift_j_3 = input('y_start: ')
            y_end = input('y_end: ')
            shift_k_3 = input('z_start: ')
            z_end = input('z_end: ')
            final_crop = filename[0:-16]+'_final_crop.nii.gz'
            shift_i_3 = int(shift_i_3)
            shift_j_3 = int(shift_j_3)
            shift_k_3 = int(shift_k_3)
            os.system('mrgrid '+f'-axis 0 {shift_i_3}:{x_end} -axis 1 {shift_j_3}:{y_end} -axis 2 {shift_k_3}:{z_end}'+' '+os.path.join(TMP, filename)+' crop '+os.path.join(TMP, final_crop))
            # shift_i_3 = -shift_i_3
            # shift_j_3 = -shift_j_3
            # shift_k_3 = -shift_k_3



            try:
                shutil.copyfile(os.path.join(TMP, final_crop),
                                os.path.join(path_crop, final_crop)
                                )
            except Exception as e:
                print('error')

    if args.Manual_Crop_coord:  
        print('\n')
        print("======================================================")
        print("=                                                    =")
        print("=               Manual Crop with txt                 =")
        print("=                                                    =")
        print("======================================================")
        print('\n')


        file_list = [x for x in os.listdir(N4_T1_path) if 'nii.gz' in x]
        for filename in file_list:
            os.remove(os.path.join(N4_T1_path, filename))
        

        file_list = [x for x in os.listdir(image_conform) if 'nii.gz' in x]
        for filename in file_list:
            os.system('cat '+image_conform+'/file_to_coords.txt grep '+os.path.join(image_conform, filename)+' '+image_conform+'/file_to_coords.txt | colrm 1 25  > '+TMP+'/coords.txt')
            f = open(os.path.join(TMP, 'coords.txt'), 'r')
            coords = f.read()
            f.close()
            size = filename[0:-7]+'_crop.txt'
            image_crop = matlab.manual_crop_brain(TMP+'/', filename, coords, size)
            f = open(os.path.join(TMP, size), 'r')
            var = f.read()
            f.close()
            print(var)
            coord = re.findall(r'\d+', var)

            final_crop = filename[0:-16]+'_final_crop.nii.gz'

            os.system('mrgrid '+var.strip('\n')+' '+os.path.join(TMP, filename)+' crop '+os.path.join(TMP, final_crop))

            try:
                shutil.copyfile(os.path.join(TMP, final_crop),
                                os.path.join(path_crop, final_crop)
                                )
            except Exception as e:
                print('error')
        



    print('\n')
    print("======================================================")
    print("=                                                    =")
    print("=                  INU Correction                    =")
    print("=                                                    =")
    print("======================================================")
    print('\n')

    file_list = [x for x in os.listdir(path_crop) if 'T1_final_crop.nii.gz' in x]
    for filename in file_list:
        N4_file = filename[0:-7]+'_N4.nii.gz'
        n4 = N4BiasFieldCorrection(output_image=os.path.join(N4_T1_path, N4_file))
        n4.inputs.dimension = 3
        n4.inputs.input_image = os.path.join(path_crop, filename)
        n4.inputs.bspline_fitting_distance = 100
        print(n4.cmdline)
        n4.run()


    print('\n')
    print("======================================================")
    print("=                                                    =")
    print("=             Average and Resampling                 =")
    print("=                                                    =")
    print("======================================================")
    print('\n')

    if not args.av_FS:

        number_images = 0
        for filenames in os.listdir(N4_T1_path):
            number_images += 1
        print("number of images: ", number_images)

        if number_images == 1:
            file_list = [x for x in os.listdir(N4_T1_path) if 'nii.gz' in x]
            for filename in file_list:
                try:
                    shutil.copyfile(os.path.join(N4_T1_path, filename),
                                    os.path.join(path_job, 'T1_preproc1.nii.gz')
                                    )
                except Exception as e:
                    print('error')
            os.system('mrgrid -voxel 0.5 '+os.path.join(path_job, 'T1_preproc1.nii.gz')+' regrid '+os.path.join(path_job, 'T1_preproc.nii.gz'))
            os.remove(os.path.join(path_job, 'T1_preproc1.nii.gz'))
        elif number_images > 1:
            os.system(f'AnatomicalAverage -s {templates_path}/NMT_05.nii.gz -o {path_job}/T1_preproc.nii.gz {N4_T1_path}/*.nii.gz && mrgrid -voxel 0.5 {path_job}/T1_preproc.nii.gz regrid {path_job}/T1_preproc.nii.gz -force')


    if args.av_FS:
        number_images = 0
        for filenames in os.listdir(N4_T1_path):
            number_images += 1
        print("number of images: ", number_images)

        if number_images == 1:
            file_list = [x for x in os.listdir(N4_T1_path) if 'nii.gz' in x]
            for filename in file_list:
                try:
                    shutil.copyfile(os.path.join(N4_T1_path, filename),
                                    os.path.join(path_job, 'T1_preproc1.nii.gz')
                                    )
                except Exception as e:
                    print('error')
            os.system('mrgrid -voxel 0.5 '+os.path.join(path_job, 'T1_preproc1.nii.gz')+' regrid '+os.path.join(path_job, 'T1_preproc.nii.gz'))
            os.remove(os.path.join(path_job, 'T1_preproc1.nii.gz'))
        elif number_images > 1:
            os.system(f'mri_motion_correct.fsl -o {path_job}/T1_preproc.nii.gz -wild *.nii.gz')
            os.system('mrgrid -voxel 0.5 '+path_job+'/T1_preproc.nii.gz regrid '+path_job+'/T1_preproc.nii.gz -force' )
            os.remove(os.path.join(path_job, 'T1_preproc.nii.gz.mri_motion_correct.fsl.log'))
            os.remove(os.path.join(path_job, 'T1_preproc.nii.gz.mri_motion_correct.fsl.log.old'))
        


    nii_in = 'T1_preproc.nii.gz'
    nii_out = 'T1_conform.nii.gz'

    data_conform = matlab.conform2(path_job+'/', nii_in, path_job+'/', nii_out)

    fslroi = ExtractROI(in_file=os.path.join(path_job, nii_in), roi_file=os.path.join(path_job, nii_in), t_min=0, t_size=1)
    fslroi.run()

    os.remove(os.path.join(path_job, 'T1_preproc.nii.gz'))

    print('\n')
    print("======================================================")
    print("=                                                    =")
    print("=                  Skull-stripping                   =")
    print("=                                                    =")
    print("======================================================")
    print('\n')

    run_Brain_Mask(os.path.join(path_job, 'T1_conform.nii.gz'), path_job+'/')
    mrconvert = mrt3.MRConvert()
    mrconvert.inputs.in_file = os.path.join(path_job, 'brain_mask.nii')
    mrconvert.inputs.out_file = os.path.join(path_job, 'mask', 'brain_mask_orig.nii.gz')
    mrconvert.run()

    os.remove(os.path.join(path_job, 'brain_mask.nii'))
    os.remove(os.path.join(path_job, 'brain.nii'))

    maths = fsl.MultiImageMaths()
    maths.inputs.in_file = os.path.join(path_job, 'mask', 'brain_mask_orig.nii.gz')
    maths.inputs.op_string = "-mul %s"
    maths.inputs.operand_files = os.path.join(path_job, 'T1_conform.nii.gz')
    maths.inputs.out_file = os.path.join(path_job, 'T1_brain.nii.gz')
    maths.run()

    maths = fsl.UnaryMaths()
    maths.inputs.in_file = os.path.join(path_job, 'mask', 'brain_mask_orig.nii.gz')
    maths.inputs.operation = 'bin'
    maths.inputs.out_file = os.path.join(path_job, 'brain_mask.nii.gz')
    maths.run()

    T1_conform = nib.load(os.path.join(path_job, 'T1_conform.nii.gz'))
    T1_conform_data = T1_conform.get_fdata()



    if not args.TMP:
        shutil.rmtree(TMP)

    

    os.system(f'flirt -in {path_job}/T1_brain.nii.gz -ref {templates_path}/NMT_template.nii -out {path_job}/{SUBID}.nii.gz -omat {path_job}/{SUBID}_acpc.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 6   -interp trilinear')
    # -schedule {FSL_DIR}/etc/flirtsch/sch3Dtrans_3dof
    os.system(f'flirt -in {path_job}/T1_conform.nii.gz -applyxfm -init {path_job}/{SUBID}_acpc.mat -out {path_job}/conform_acpc.nii.gz -paddingsize 0.0 -interp trilinear -ref {path_job}/{SUBID}.nii.gz')
    
    os.system(f'gzip -d {path_job}/conform_acpc.nii.gz')

    # file_list = [x for x in os.listdir(T1_image_path) if 'nii.gz' in x]
    # for filename in file_list:
    try:
        shutil.copyfile(T1_image_path,
                        os.path.join(path_job, 'T1.nii.gz')
                        )
    except Exception as e:
        print('error')



    os.system(f'gzip -d {path_job}/T1.nii.gz')

    if args.Sphinx:
        cvt = MRIConvert()
        cvt.inputs.in_file = os.path.join(path_job, 'T1.nii')
        cvt.inputs.out_file = os.path.join(path_job, 'T1.nii')
        cvt.inputs.sphinx = True
        cvt.run()



    os.system(f'gzip -d {path_job}/simu/CT.nii.gz')
    os.system(f'gzip -d {path_job}/T1_conform.nii.gz')
    print('> Coords of AC(world space)')
    os.system(f'fsleyes {path_job}/{SUBID}.nii.gz')
    
    acpc_i = input('i: ')
    acpc_j = input('j: ')
    acpc_k = input('k: ')
    
   
    if args.asym:
        os.system(f'bash align2NMT.sh -id {SUBID} -o {OUTPUT_DIR} -t ./NMT_v2.0_asym/NMT_v2.0_asym_05mm')
    else:
        os.system(f'bash align2NMT.sh -id {SUBID} -o {OUTPUT_DIR} -t ./NMT_v2.0_sym/NMT_v2.0_sym_05mm')

    
    file_list = [x for x in os.listdir(path_job) if 'nii.gz' in x]
    for filename in file_list:
        os.system(f'gzip -d {path_job}/{filename}')

    # os.system(f'gzip -d {path_job}/{SUBID}.nii.gz')
    # os.system(f'gzip -d {path_job}/SEG_in_{SUBID}.nii.gz')
    
    
    matlab.acpc_job(f'{path_job}/{SUBID}.nii', acpc_i, acpc_j, acpc_k)
    matlab.acpc_job(f'{path_job}/conform_acpc.nii', acpc_i, acpc_j, acpc_k)

    # 由于MRI图像的FOV不足以概括整个头部，即使之前的分割程序已经拓宽成了256*256*256的矩阵。因此将CT配准到MRI图像之前，需要对T1_conform图像进行padding
    # matlab.padding(f'{path_job}', f'SEG_in_{SUBID}.nii')
    matlab.padding(f'{path_job}')
    file_list = [x for x in os.listdir(path_job) if f'_in_{SUBID}' in x]
    for filename in file_list:
        matlab.padding_all(f'{path_job}', filename, 'simu/'+filename[0:-4]+'_padding.nii')


    # 将原始T1图像的坐标配准到conform_acpc
    matlab.t1_job(f'{path_job}')
    # 将CT图像配准到T1
    matlab.ct_job(f'{path_job}')
    # 基于T1_padding对CT图像进行reslice
    matlab.reslice_job(f'{path_job}')
    matlab.exit()

    t1 = time.time()
    print('\n')
    print("TOTAL RUNNING TIMES: ", t1-t0, "seconds")
