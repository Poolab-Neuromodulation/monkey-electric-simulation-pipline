# Prerequisite

It runs only on Linux system

## Tools

You previously need the installation of:

 - **FREESURFER** v6.0.0.1. https://surfer.nmr.mgh.harvard.edu/pub/dist/freesurfer/  
 - **FSL**. https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation 
 - **ANTS**. http://stnava.github.io/ANTs/
 - **MATLAB**. https://la.mathworks.com/products/get-matlab.html?s_tid=gn_getml 

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

5. Set the electrode position and current level you want to stimulate, and save it as a json file.

   ```json
   {
      "Fpz" : 1, 
      "T1" : -1
   }
   ```

6. run the pipeline:

   ```shell
   python run.py -t1 {T1_image_path} -ct {CT_image_path} -id {SUBID} -o {OUTPUT_DIR} -e {json_path}
   ```
