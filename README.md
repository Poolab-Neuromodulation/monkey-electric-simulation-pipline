# Prerequisite

It runs only on Linux system

## Tools

You previously need the installation of:

 - **FREESURFER** v6.0.0.1. https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/  
 - **FSL**. https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation 
 - **ANTS**. http://stnava.github.io/ANTs/
 - **AFNI**. https://afni.nimh.nih.gov/
 - **MATLAB**. https://la.mathworks.com/products/get-matlab.html?s_tid=gn_getml   
 -- MATLAB Toolbox: spm12, iso2mesh, NIfTI_Tools_matlab

Make sure the above tools are added to the environment variables.

## python environment

* **anaconda** (LInux). https://www.anaconda.com/

```shell
conda env create -f environment.yml
```
# How to use

1. Convert dcm files of raw CT and T1  to nii files

2. If CT and T1 images are not in sphinx position, do reorientation:

   ```shell
   mri_convert CT.nii --sphinx CT.nii
   mri_convert T1.nii --sphinx T1.nii
   ```

3. ```shell
   gzip CT.nii # nii => nii.gz
   gzip T1.nii
   ```

4. alignment to anterior and posterior commissure for CT and T1


5. run the pipeline:

   ```shell
   python run.py -t1 {T1_image_path} -ct {CT_image_path} -id {SUBID} -o {OUTPUT_DIR}
   ```

# Environment Package
链接：https://pan.baidu.com/s/186mmE3rDBHBEIrXlSzrzeg?pwd=cr91
提取码：cr91

1. install apt package

   ```shell
   sudo apt-get update
   sudo apt-get -y install itksnap
   sudo apt-get -y install bc binutils libgomp1 perl psmisc sudo tar tcsh unzip uuid-dev vim-common libjpeg62-dev
   sudo bash OS_notes.linux_ubuntu_20_64_a_admin.txt 2>&1 | tee o.ubuntu_20_a.txt
   sudo apt install cmake
   sudo apt -y install git
   ```
   
2. install anaconda
   ```shell
   echo -e "\nyes\n\nyes\n" | bash ./Anaconda3-2021.11-Linux-x86_64.sh 
   source ~/.bashrc
   ```

3. install all packages except Matlab
   ```shell
   bash ./environment_install.sh 
   ```
