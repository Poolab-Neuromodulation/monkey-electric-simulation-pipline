
function reslice_job(dn)
    matlabbatch = reslice(dn);
    spm('defaults', 'PET');
    spm_jobman('run', matlabbatch);
end

function matlabbatch = reslice(dn)

	source_image = fullfile(dn, 'simu/CT.nii,1')
	ref_image = fullfile(dn, 'T1_padding.nii,1')
	
	matlabbatch{1}.spm.spatial.coreg.write.ref = {ref_image};
	matlabbatch{1}.spm.spatial.coreg.write.source = {source_image};
	matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
	matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
	matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
	matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = '_reslice';
	
end    
