function smoothAndWarpContrasts(rootData,sub,glm,fwhm)
% smooth the contrasts using SPM.  
% then use the transformations calculated by FSL (during preprocessing) to
% warp the contrasts to standard space. If you're using SPM also for
% preprocessing this will look different. 

% fwhm: value in mm of fwhm of smoothing kernel

spmPath = '/home/fs0/abaram/scratch/MATLAB/spm12';
addpath(genpath(spmPath));

%% smooth contrasts 

spmDir = fullfile(rootData,'spmDirs',glm,sub);
regDir = fullfile(rootData,'reg',sub); % this is the /reg dir from the FSL preprocessing folder, with transformation computed by FLIRT and FNIRT
smthConDir = fullfile(spmDir,['smth' num2str(fwhm)]);
cd(spmDir)
clear matlabbatch
load SPM;

% find all contrast images of the current subject
conImages  = spm_select('List',['con*nii']);
for iCon = 1:size(conImages,1)
    matlabbatch{1}.spm.spatial.smooth.data{iCon,1} = fullfile(spmDir,conImages(iCon,:));
end
% specify smoothing kernel fwhm
matlabbatch{1}.spm.spatial.smooth.fwhm = [fwhm fwhm fwhm];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = ['smth' num2str(fwhm) '_'];

% smooth contrast images, place in glm folder
spm_jobman('run',matlabbatch);
movefile([spmDir,'/smth' num2str(fwhm) '_*'],smthConDir);

%% warp contrasts

mkdir(fullfile(smthConDir,'warp'));
cd(fullfile(smthConDir));
imagesToWarp = spm_select('List',['smth' num2str(fwhm) '_con*nii']); %  all smoothed contrasts
warpedImagesFileNames = [imagesToWarp(:,1:end-4),repmat('_MNI',size(imagesToWarp,1),1)]; % add _MNI to warped file names
%% use FSL to warp contrast to standard space
for iCon = 1:size(imagesToWarp,1)
   system(['applywarp -i ' imagesToWarp(iCon,:) ' -o ' fullfile(smthConDir,'warp',warpedImagesFileNames(iCon,:)) '.nii.gz -r ' regDir '/standard.nii.gz -w ' regDir '/highres2standard_warp.nii.gz']);
end
