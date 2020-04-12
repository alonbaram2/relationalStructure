function searchlightDefinitionVolumeWrapper(rootData,sub,roi)
% Define searchlights within a volumetric RPOI. This is done in MNI space,
% as otherwise it's difficult later to match between grey matter of
% different subjects. 

% This was used on the ROI "lAccumbens20": The left Accumbens label from the
% Harvard-Oxford subcortical structural atlas, thresholded at a probability
% of 20 (i.e. any voxels that are part of the left nucleus accumbens with a
% probabiliyt of more than 0.2). 

anatomicalDir = fullfile(rootData,'anatomical','mni152Masks');

% mask to find searchlights in. Only the voxels corresponding to the
% cortical surface within this mask will be used. 
Vmask = spm_vol(fullfile(anatomicalDir,[roi '.nii'])); % mask
Vmask.data = spm_read_vols(Vmask);
Vmask.sphere=[10,100];

% File name of searchlight matlab structure
fname = ['SL_vol_' roi '_v100.mat']; % 100 voxels in each searchlight
slDir = fullfile(rootData,'searchlights',sub);
if ~exist(slDir,'dir') % create a folder for 1st-level results
    mkdir(slDir);
end

L = rsa.defineSearchlight({Vmask},Vmask);
save(fullfile(slDir,fname),'-struct','L'); 