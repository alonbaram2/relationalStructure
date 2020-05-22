rootData    = '/home/fs0/abaram/scratch/NN2020';

rootScripts = '/home/fs0/abaram/scripts/NN2020';
addpath(genpath(rootScripts));
spmPath = '/home/fs0/abaram/scratch/MATLAB/spm12';
addpath(genpath(spmPath));

sub='sub-00q';
glm='GLM1';
roi='surf';
fwhm = 5;
iBlock = 999;
% assuming already ran getEvs.m for current subject

% % run SPM first-level GLM
% runSpmGlmEst(rootData,sub,'confounds')
% runRSA(rootData,sub,'GLM2','xRun','surf',fwhm)


subjects = {'sub-01','sub-02','sub-03','sub-04','sub-05','sub-06','sub-07','sub-08','sub-09',...
    'sub-10','sub-11','sub-12','sub-13','sub-14','sub-15','sub-16','sub-17','sub-18','sub-19',...
    'sub-20','sub-21','sub-22','sub-23','sub-24','sub-25','sub-26','sub-27','sub-28'};

runSpmGlmEst(rootData,sub,'confounds')
confound = getConvolvedConfoundsMat(rootData,subjects,str2num(sub(end-1:end)));
confound.name = 'allConfounds';
runRSA_confounds(rootData,sub,'GLM2','xRun','surf',fwhm,confound)

% % run SPM first-level GLM
% runSpmGlmEst_MNI(rootData,sub,'GLM2')

% runRSA(rootData,sub,glm,'xRun',roi,5)



% runRSA_withinCluster(rootData,sub,'GLM2','corrWithinBlock',roi,5,iBlock)
