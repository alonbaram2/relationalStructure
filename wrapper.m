clear

rootScripts = '/home/fs0/abaram/scripts/NN2020';
addpath(genpath(rootScripts));
spmPath = '/home/fs0/abaram/scratch/MATLAB/spm12';
addpath(spmPath);

rootData    = '/home/fs0/abaram/scratch/NN2020_BIDS/derivatives';

subjects = {'sub-01','sub-02','sub-03','sub-04','sub-05','sub-06','sub-07',...
    'sub-11','sub-12','sub-13','sub-14','sub-15','sub-16','sub-17','sub-18','sub-19',...
    'sub-20','sub-21','sub-22','sub-23','sub-24','sub-25','sub-26','sub-27','sub-28'};

nSub = length(subjects);
nRun  = 2; % two independent runs, each with the 4 conditions
nCond = 4; % the 4 experimental conditions: 2 correlation typs, 2 stimuli sets.
nTrialsPerBlock = 30;

fwhm = 5; % smoothing fwhm

%% Univariate results (Fig 1e and S3): Chosen action value, with STRUCT and NAIVE models competing
% ###  GLM1 - first level
for iSub=1:length(subjects)
    % get regressors
    getEvs(rootData,subjects,iSub,'GLM1')
    
    % run SPM first-level GLM
    runSpmGlmEst(rootData,subjects{iSub},'GLM1')
    
    % calculate contrasts
    runSpmContrast(rootData,subjects{iSub},'GLM1')
    
    % warp contrasts to standard space and smooth
    smoothAndWarpContrasts(rootData,subjects{iSub},'GLM1')   
end

% second level GLM - Fig 1e
univariateGroupLevel(rootData,subjects,'GLM1')

%% For all RSA analyses: define the searchlight on the surface. 
% This code assumes you have already ran the Freesurfer surface reconstruciton
% using the recon-all command: https://surfer.nmr.mgh.harvard.edu/fswiki/recon-all
% This code is adapted from the RSA toolbox, and requires it to be in your
% path: https://github.com/rsagroup/rsatoolbox
for iSub=1:length(subjects)    
    searchlightDefinitionSurfaceWrapper(rootData,subjects{iSub})
end

%% Representation of relational structure between stimuli (Fig 2)
% ### GLM2 - first level

for iSub=1:length(subjects)
    % get regressors
    getEvs(rootData,subjects,iSub,'GLM2')
    
    % run SPM first-level GLM
    runSpmGlmEst(rootData,subjects{iSub},'GLM2')
end         

% run the RSA analysis on the cortical surface. 
for iSub=1:length(subjects)
    % get the extended RDM (i.e. before crossing the blocks)
    runRSA(rootData,subjects{iSub},'GLM2','xRun','surf',fwhm)
end

% # group analysis 
% first stack all single subjects contrasts into one file
stackSubjectsContrasts(rootData,subjects,'GLM2','xRun','surf');

% Run permutation tests - Fig 2d and S4
    % run permutation tests in PALM on the surface
    % The permutation tests also correct for multiple comparisons. 
    % Note that  this takes a long time to run.
    clusterThresh = '3.1'; % cluster forming threshold. 
    nPerm       = '10000';% number of permutations
    maskName = 'CORTEX';
    % run PALM
    groupLevelPermTests(rootData,'GLM2','xRun','surf',clusterThresh,nPerm,maskName)
    
% # plot RDMs at the peaks: entorhinal for relational structure and 
%   LOC for visual identity  (both in right hemisphere)
%   First get vertex of peak activation, then plot average RDM across
%   participants in this vertex. 

% visualise effects at their peak
% relational structure - Fig 2b and c, top. Get the max vertex from the group level map with
% all 28 subjects, saved in masksAndRois. This the 28-subject result equivalent to
% the 25-subjects files in fullfile(rootData,'RDMs','GLM2','groupStats','perm','xRun','surf','CORTEX')
relationalStructMap = load_mgh(fullfile(rootData,'masksAndRois','fsaverage','relationalStructure_smth5_rh_nPerm10000_clstrTh3p1_dpv_tstat_uncp_c1_all28subjects.mgz'));
[~,indMax] = max(relationalStructMap);
roiStr = [num2str(indMax - 1),'rh']; % -1 is for conversion from Matlab to Freesurfer coords. 
plotDataRdm(rootData,subjects,'GLM2','xRun','relationalStructure',roiStr,fwhm)
plotGardnerAltman(rootData,subjects,'GLM2','xRun','surf','relationalStructure',roiStr,fwhm)

% # visual identity - Fig 2b and c, bottom
visualIdentityMap = load_mgh(fullfile(rootData,'masksAndRois','fsaverage','visualIdentity_smth5_rh_nPerm10000_clstrTh3p1_dpv_tstat_uncp_c1_all28subjects.mgz'));
[~,indMax] = max(visualIdentityMap);
roiStr = [num2str(indMax - 1),'rh']; % -1 is for conversion from Matlab to Freesurfer coords. 
plotDataRdm(rootData,subjects,'GLM2','xRun','visualIdentity',roiStr,fwhm)  
plotGardnerAltman(rootData,subjects,'GLM2','xRun','surf','visualIdentity',roiStr,fwhm)

%% Prediction error x structure interaction (Fig 3)
% GLM3 - first level. 
for iSub=1:length(subjects)
    % get regressors
    getEvs(rootData,subjects,iSub,'GLM3')
    
    % run SPM first-level GLM
    runSpmGlmEst(rootData,subjects{iSub},'GLM3')
    
    % calculate contrast
    runSpmContrast(rootData,subjects{iSub},'GLM3')
    
    % warp contrast to standard space and smooth
    smoothAndWarpContrasts(rootData,subjects{iSub},'GLM3')    
end
% second level GLM - the peaks of this contrast (ventral striatum and
% vmPFC) will serve as ROIs for the PE x structure interaction effect.
% Fig 3b and c (insets). 
univariateGroupLevel(rootData,subjects,'GLM3')
projectTMapsToSurface(rootData,'GLM3')

% run the RSA analysis on the cortical surface. 
% To perform statistical tests we will perform an ROI analysis on the peaks
% of the univatiate STRUCT prediction error effect from the same GLM (GLM3). 
for iSub=1:length(subjects)
    runRSA(rootData,subjects{iSub},'GLM3','xRun','surf',fwhm)
end

% # group analysis 
% first stack all single subjects contrasts into one file
stackSubjectsContrasts(rootData,subjects,'GLM3','xRun','surf');
% Run permutation tests in PALM.
% The permutation tests also correct for multiple comparisons. %
% However this takes a long time to run.
clusterThresh = 'None'; % cluster forming threshold. For this GLM we are 
                        % only looking at ROIs defined by the peaks of the
                        % univarite prediction error effect. 
nPerm       = '10000';% number of permutations
% mask to run the permutation tests in and perform multiple comparisons in.
maskName = 'CORTEX';
% run PALM
groupLevelPermTests(rootData,'GLM3','xRun','surf',clusterThresh,nPerm,maskName)

% get ROI for vmPFC effect from univariate PE effect: note I'm using the
% T-map from all 28 subjects (and not just the 25 subjects who agreed to
% share their data). This is provided in rootData/masksAndRois/fsaverage.
% this is the 28 subjects equivalent of
% rootData/spmDirs/GLM3/2ndLevel_all/crPe_STRUCT/spmT_GLM3_crPe_STRUCT_0001_lh.mgh

peUnivariateSurfaceEffect_lh = load_mgh(fullfile(rootData,'masksAndRois','fsaverage','spmT_GLM3_crPe_STRUCT_0001_lh_all28subjects.mgh'));
peUnivariateSurfaceEffect_rh = load_mgh(fullfile(rootData,'masksAndRois','fsaverage','spmT_GLM3_crPe_STRUCT_0001_rh_all28subjects.mgh'));
[max_lh,roiVertex_lh] = max(peUnivariateSurfaceEffect_lh);
[max_rh,roiVertex_rh] = max(peUnivariateSurfaceEffect_rh);
roiStr_lh = [num2str(roiVertex_lh - 1) 'lh']; % -1 is for conversion to 0-based Freesurfer indeces
roiStr_rh = [num2str(roiVertex_rh - 1) 'rh']; % -1 is for conversion to 0-based Freesurfer indeces
% plot subject-average data RDM in the vmPFC peaks of univariate effect: Fig
% 3b,c and S5
plotGardnerAltman(rootData,subjects,'GLM3','xRun','surf','peXstructure',roiStr_lh,fwhm)
plotGardnerAltman(rootData,subjects,'GLM3','xRun','surf','peXstructure',roiStr_rh,fwhm)
plotDataRdm(rootData,subjects,'GLM3','xRun','peXstructure',roiStr_lh,fwhm)

% # run a parallel pipeline in MNI space rather than subject space, in
% order to run RSA on volumetric searchlights (specifically in the ventral
% Striatum - whole brain peak of the PE univariate effect). 

for iSub=1:length(subjects)
    % run SPM first-level GLM in MNI space
    runSpmGlmEst_MNI(rootData,subjects{iSub},'GLM3')
    % define volumetric searchlights only within the left nucleus
    % accumbens, as that's where the univariate "correctness prediction error"
    % effect peaks. This will be used as an ROI for the PE x structure
    % interaction effect, together with the vmPFC univariate peal 
    searchlightDefinitionVolumeWrapper(rootData,subjects{iSub},'lAccumbens20') % 
    runRSA(rootData,subjects{iSub},'GLM3','xRun','lAccumbens20',fwhm) 
end

% group stats
% group stats on volumetric data in ventral striatum 
% first stack all single subjects contrasts into one file
stackSubjectsContrasts(rootData,subjects,'GLM3','xRun','lAccumbens20');
clusterThresh = 'None'; % searching only in a peak voxel (from the orthogonal univariate contrast).
nPerm       = '10000';% number of permutations
maskFile          = fullfile(rootData,'masksAndRois','mni','lAccumbens20.nii');
groupLevelPermTests(rootData,'GLM3','xRun','lAccumbens20',clusterThresh,nPerm,maskFile)

% Fig 3c: peak of univariate PE effect in all the brain is ventral striatum
peXStructureUnivariateSurfaceEffect_vol = read_avw(fullfile(rootData,'masksAndRois','fsaverage','spmT_GLM3_crPe_STRUCT_0001_all28subjects.nii'));
[~,roiVox] = max(peXStructureUnivariateSurfaceEffect_vol(:));
[x,y,z] = ind2sub(size(peXStructureUnivariateSurfaceEffect_vol),roiVox); 
% plot interaction effect at peak of univariate effect
plotGardnerAltman(rootData,subjects,'GLM3','xRun','lAccumbens20','peXstructure',[x,y,z] - 1,fwhm) % -1 is for conversion to 0-based FSL indeces

%% Representation of a behaviourally relevant latent variable. (Fig 4)
% For this analysis we also use GLM2, only here we will use all 5 events
% ([A], [B], [C] during both stimulus and outcome presentation

% get single subject contrasts. note that the distType is now "corr" - we
% are looking at within-block effects, measures with simple correlation
% distances (with spatial whitening). 
for iSub=1:length(subjects)
    runRSA(rootData,subjects{iSub},'GLM2','corrWithinBlock','surf',fwhm)
end

% # group analysis - Fig 4a right
% first stack all single subjects contrasts into one file
stackSubjectsContrasts(rootData,subjects,'GLM2','corrWithinBlock','surf');
% Run permutation tests
% run on the surface permutation tests in PALM.
clusterThresh = '2.3'; % cluster forming threshold.
nPerm       = '10000';% number of permutations
maskName = 'PALS_B12_OFC'; % change this to 'CORTEX' if you want to look in all cortex (e.g. for the pPHCg effect).
% run PALM
groupLevelPermTests(rootData,'GLM2','corrWithinBlock','surf',clusterThresh,nPerm,maskName)

% get index of max of latent variable effect. Once again get the index from
% the group map of all 28 subjects, provided at
% rootData/masksAndRois/fsaverage/
LatentVariableMap = load_mgh(fullfile(rootData,'masksAndRois','fsaverage',...
                        'latentVar_all_smth5_lh_nPerm10000_clstrTh2p3_dpv_tstat_c1_all28subjects.mgz'));
[~,indMax] = max(LatentVariableMap);
roiStr = [num2str(indMax - 1),'lh']; % -1 is for conversion from Matlab to Freesurfer coords. 
% Fig 4a left and 4b
plotGardnerAltman(rootData,subjects,'GLM2','corrWithinBlock','surf','latentVar_all',roiStr,fwhm)
plotGardnerAltman(rootData,subjects,'GLM2','corrWithinBlock','surf','latentVar_S-S',roiStr,fwhm)
plotGardnerAltman(rootData,subjects,'GLM2','corrWithinBlock','surf','latentVar_O-O',roiStr,fwhm)
plotGardnerAltman(rootData,subjects,'GLM2','corrWithinBlock','surf','latentVar_S-O',roiStr,fwhm)

% plot RDMs at the peak of OFC effect  - Fig S6b. 
plotDataRdm(rootData,subjects,'GLM2','corrWithinBlock','latentVar_S-S',roiStr,fwhm)
plotDataRdm(rootData,subjects,'GLM2','corrWithinBlock','latentVar_O-O',roiStr,fwhm)
plotDataRdm(rootData,subjects,'GLM2','corrWithinBlock','latentVar_S-O',roiStr,fwhm)


