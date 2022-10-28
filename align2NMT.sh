for arg in "$@"
do
  case "$arg" in
  -id)
   SUB_ID=$2
   shift;shift
  ;;
  -o)
   Output=$2
   shift;shift
  ;;
  -t)
   NMT_Path=$2
   shift;shift
  ;;
   esac
done


CHARM_Path=${NMT_Path}/supplemental_CHARM
SARM_Path=${NMT_Path}/supplemental_SARM

@animal_warper -input ${Output}/${SUB_ID}/${SUB_ID}.nii.gz -base ${NMT_Path}/NMT_*_SS.nii.gz -base_abbrev NMT2 -atlas_followers ${CHARM_Path}/CHARM_6_*.nii.gz ${CHARM_Path}/CHARM_5_*.nii.gz ${CHARM_Path}/CHARM_4_*.nii.gz ${CHARM_Path}/CHARM_3_*.nii.gz ${CHARM_Path}/CHARM_2_*.nii.gz ${CHARM_Path}/CHARM_1_*.nii.gz ${SARM_Path}/SARM_6_*.nii.gz ${SARM_Path}/SARM_5_*.nii.gz ${SARM_Path}/SARM_4_*.nii.gz ${SARM_Path}/SARM_3_*.nii.gz ${SARM_Path}/SARM_2_*.nii.gz ${SARM_Path}/SARM_1_*.nii.gz -atlas_abbrevs CHARM_6 CHARM_5 CHARM_4 CHARM_3 CHARM_2 CHARM_1 SARM_6 SARM_5 SARM_4 SARM_3 SARM_2 SARM_1 -seg_followers ${NMT_Path}/NMT_*_segmentation.nii.gz -seg_abbrevs SEG -outdir ${Output}/${SUB_ID}/ -ok_to_exist 

# -skullstrip ${NMT_Path}/NMT_*_brainmask.nii.gz

