rootData    = '/home/fs0/abaram/scratch/NN2020';

rootScripts = '/home/fs0/abaram/scripts/NN2020';
addpath(genpath(rootScripts));
spmPath = '/home/fs0/abaram/scratch/MATLAB/spm12';
addpath(genpath(spmPath));

sub='sub-01';
glm='GLM1';
roi='surf';
fwhm = 5;
iBlock = 999;
% assuming already ran getEvs.m for current subject

% % run SPM first-level GLM
% runSpmGlmEst(rootData,sub,'GLM2')
% runRSA(rootData,sub,'GLM2','xRun','surf',fwhm)


% 
% % calculate contrasts
% runSpmContrast(rootData,sub,'GLM1')
% 
% % 
% smoothAndWarpContrasts(rootData,sub,'GLM1')
% 
% %
% runRSA(rootData,sub,'GLM2','xRun','surf')


% % run SPM first-level GLM
% runSpmGlmEst_MNI(rootData,sub,'GLM2')

% runRSA(rootData,sub,glm,'xRun',roi,5)


runRSA_withinCluster(rootData,sub,'GLM2','corrWithinBlock',roi,5,iBlock)
