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

1. ac pc校正
2. 

```shell
# 1
python preprocess.py -t1 {T1_image_path} -ct {CT_image_path} -id {SUBID} -o {OUTPUT_DIR}

# 2
python ct_segmentation.py -id {SUBID} -o {OUTPUT_DIR}

# 3
python simu_monkey.py -id {SUBID} -o {OUTPUT_DIR} -e {json_path}

# 4 all
python run.py -t1 {T1_image_path} -ct {CT_image_path} -id {SUBID} -o {OUTPUT_DIR} -e {json_path}
```

