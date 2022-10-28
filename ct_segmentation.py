from cv2 import convertFp16
from nbformat import read
import numpy as np
import nibabel as nib
from sklearn.cluster import DBSCAN
import cv2
import argparse, sys, os

def load_arguments():
    parser = argparse.ArgumentParser(description="Remove bed plate for CT image")
    parser.add_argument("-o", "--output", help="Directory where the output file is to be saved")
    parser.add_argument("-id", "--ID", help="SUBJECT ID", type=str)
    args = parser.parse_args()
    return args


def normalize(image):
    mean = np.mean(image)
    var = np.mean(np.square(image-mean))
    image = (image - mean)/np.sqrt(var)
    return image
    
def read_nii(nii_path):
    nii_obj = nib.load(nii_path)
    nii_data = nii_obj.get_fdata()
    affine = nii_obj.affine
    hdr = nii_obj.header

    return nii_data, affine, hdr

if __name__ == "__main__":
    args = load_arguments()

    if not args.output:
        print("Should specify output. See -h for instructions.")
        sys.exit(1)
    
    if not args.ID:
        print("Should specify subject ID. See -h for instructions.")
        sys.exit(1)

    # PATH
    OUTPUT_DIR = args.output # out_path
    SUBID = args.ID

    print("removing bed from CT image\n")

    ct_path = os.path.join(OUTPUT_DIR, SUBID, 'simu')
    os.system(f'gzip {ct_path}/CT_reslice.nii')

    # ct_image = f'{ct_path}/CT_reslice.nii.gz'
    # ct_obj = nib.load(ct_image)
    # ct_data = ct_obj.get_fdata()
    # affine = ct_obj.affine
    # hdr = ct_obj.header
    ct_data, affine, hdr = read_nii(f'{ct_path}/CT_reslice.nii.gz')

    # 让ct图像中padding部分设为-1000
    img = np.where(np.isnan(ct_data), -1000, ct_data)
    img = np.where(img == 0, -1000, img)
    img = np.where(img < -1000, -1000, img)

    # 标准化
    std_img = normalize(img)
    new_nii = nib.Nifti1Image(std_img, affine, hdr) 
    nib.save(new_nii, f'{ct_path}/norm_ct.nii.gz')


    for i in range(len(std_img[0][0])):
        img = std_img[:, :, i]
    
        ret, binary = cv2.threshold(np.uint8(img), 0, 255, cv2.THRESH_BINARY) 
        # binary = cv2.adaptiveThreshold(np.uint8(img), 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY,11, 4)

        std_img[:, :, i] = binary

    train = (std_img == 255)
    x = np.unique(np.where(std_img==255)[0])

    for i in x:
        coord = []
        for j in range(len(std_img[i])):
            for k in range(len(std_img[i][j])):
                if train[i][j][k] == True:
                    coord.append([j, k])
        coord = np.asarray(coord)  

        y_pred = DBSCAN(eps=3).fit(coord).labels_

        label_num = np.unique(y_pred)

        bed_label = y_pred[0]
        bed_pos = np.where(y_pred == bed_label)[0][0]
        # print(bed_pos)

        for n in range(len(label_num)):
            """ 遍历DBSCAN聚类所有类别, 不管分几类，床板的位置一定在最下面，z坐标较小的那个类别就是床板的类别 """
            pos = np.where(y_pred == label_num[n])[0][0]

            if coord[pos][1] < coord[bed_pos][1]:
                bed_pos = pos

        bed_label = y_pred[bed_pos]
        for m in np.where(y_pred == bed_label)[0]:
            std_img[i][coord[m][0]][coord[m][1]] = 127.0

    os.system(f'gzip -d {ct_path}/CT_reslice.nii')
    bed_img = np.where(std_img == 127, 1, 0)
    new_nii = nib.Nifti1Image(std_img, affine, hdr) 
    nib.save(new_nii, f'{ct_path}/ct_seg.nii.gz')

    new_nii = nib.Nifti1Image(bed_img, affine, hdr) 
    nib.save(new_nii, f'{ct_path}/bed.nii.gz')

    os.system(f'fslmaths {ct_path}/ct_seg.nii.gz -bin {ct_path}/ct_seg.nii.gz')
    os.system(f'fslmaths {ct_path}/ct_seg.nii.gz -sub {ct_path}/bed.nii.gz -bin {ct_path}/skin.nii.gz')


    skin_data, affine, hdr = read_nii(f'{ct_path}/skin.nii.gz')
    ret, th = cv2.threshold(np.uint8(skin_data), 0, 255, cv2.THRESH_BINARY)

    im_floodfill = th.copy()

    for i in range(len(skin_data)):
        h, w = skin_data[:, i, :].shape[:2]
        mask = np.zeros((h+2, w+2), np.uint8)

        cv2.floodFill(im_floodfill[:, i, :], mask, (0,0), 255)
        im_floodfill[:, i, :] = cv2.bitwise_not(im_floodfill[:, i, :])

    new_nii = nib.Nifti1Image(im_floodfill, affine, hdr) 
    nib.save(new_nii, f'{ct_path}/air.nii.gz')

    os.system(f'fslmaths {ct_path}/air.nii.gz -bin {ct_path}/air.nii.gz')
    # os.system(f'fslmaths {ct_path}/skin.nii.gz -add {ct_path}/air.nii.gz -bin -dilM {ct_path}/skin.nii.gz')
    # os.system(f'fslmaths {ct_path}/skin.nii.gz -sub {ct_path}/air.nii.gz -bin -dilM {ct_path}/skin.nii.gz')

    norm_data, affine, hdr = read_nii(f'{ct_path}/norm_ct.nii.gz')
    norm_img = np.uint8(np.where(norm_data < 0, 0, norm_data))
    for i in range(len(norm_img[0][0])):
        img = norm_img[:, :, i]

        ret, binary = cv2.threshold(np.uint8(img), 2, 1, cv2.THRESH_BINARY) 
        # binary = cv2.adaptiveThreshold(np.uint8(img), 1, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 7, 5)

        norm_img[:, :, i] = binary
    
    new_nii = nib.Nifti1Image(norm_img, affine, hdr) 
    nib.save(new_nii, f'{ct_path}/bone.nii.gz')

    os.system(f'fslmaths {ct_path}/bone.nii.gz -sub {ct_path}/bed.nii.gz -bin {ct_path}/bone.nii.gz')

    parcel_img, affin, hdr = read_nii(f'{ct_path}/SEG_in_{SUBID}_padding.nii') 
    bones_img, affin, hdr = read_nii(f'{ct_path}/bone.nii.gz') 
    air_img, affin, hdr = read_nii(f'{ct_path}/air.nii.gz') 
    skin_img, affin, hdr = read_nii(f'{ct_path}/skin.nii.gz')

    label_tissue = np.zeros(skin_img.shape)

    # wm
    label_tissue += np.where(parcel_img == 4, 1, 0)
    label_tissue += np.where(parcel_img == 5, 1, 0)
    label_tissue = np.where(label_tissue > 1, 1, label_tissue)
    # gm
    label_tissue += np.where(parcel_img == 2, 2, 0) 
    label_tissue += np.where(parcel_img == 3, 2, 0)
    label_tissue = np.where(label_tissue > 2, 2, label_tissue)
    # csf
    label_tissue += np.where(parcel_img == 1, 3, 0)
    label_tissue = np.where(label_tissue > 3, 3, label_tissue)
    # bones
    label_tissue += np.where(bones_img == 1, 4, 0)
    label_tissue = np.where(label_tissue > 4, 4, label_tissue)
    # skin
    label_tissue += np.where(skin_img == 1, 5, 0)
    label_tissue = np.where(label_tissue > 5, label_tissue - 5, label_tissue)
    # air
    label_tissue += np.where(air_img == 1, 6, 0)
    label_tissue = np.where(label_tissue > 6, 5, label_tissue)


    new_nii = nib.Nifti1Image(label_tissue, affine, hdr) 
    nib.save(new_nii, f'{ct_path}/label_tissue.nii.gz')

    os.remove(f'{ct_path}/norm_ct.nii.gz')
    os.remove(f'{ct_path}/ct_seg.nii.gz')
    os.system(f'gzip -d {ct_path}/label_tissue.nii.gz')