function write_SL_mask(vox,sub,sl)

% vox: 3x1 center of SL voxel indeces as appear in FSLeyes. must be in
% cortical grey matter.

rootData = '/vols/Scratch/abaram/NN2020';
if strcmp(sl,'surf')
    SLfname = fullfile(rootData,'searchlights',sub,['SL_' sl '_v100.mat']);
    V = spm_vol([rootData '/anatomical/' sub '/structural_brain_mask_2mm.nii']);
else 
    SLfname = fullfile(rootData,'searchlights',sub,['SL_vol_' sl '_v100.mat']);
    V = spm_vol([rootData '/anatomical/mni152Masks/brain.nii']);
    
end
load(SLfname,'LI','voxel');
linIndCenter = sub2ind(V.dim,vox(1)+1,vox(2)+1,vox(3)+1); % add 1 for conversion from FSLeyes to Matlab

if ~any(linIndCenter==voxel)
    error('center voxel is not in mask')
end
mask = zeros(V.dim);
mask(LI{linIndCenter==voxel}) = 1;
save4Dnii(mask,V,fullfile(rootData,'searchlights',sub, ['vox' num2str(vox(1)) num2str(vox(2)) num2str(vox(3)) '_SLmask_' sl '.nii']))
