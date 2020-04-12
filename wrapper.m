clear

rootScripts = '/home/fs0/abaram/scripts/NN2020';
addpath(genpath(rootScripts));
spmPath = '/home/fs0/abaram/scratch/MATLAB/spm12';
addpath(genpath(spmPath));

groupStatsPermTestsFlag = false; % run permutation tests at the group level for RSA- very long

rootData    = '/home/fs0/abaram/scratch/NN2020';
subjects = {'sub-01','sub-02','sub-03','sub-04','sub-05','sub-06','sub-07','sub-08','sub-09',...
    'sub-10','sub-11','sub-12','sub-13','sub-14','sub-15','sub-16','sub-17','sub-18','sub-19',...
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
stackSubjectsContrasts(rootData,'GLM2','xRun','surf');

% Run permutation tests - Fig 2d and S4
    % run permutation tests in PALM on the surface
    % The permutation tests also correct for multiple comparisons. 
    % Note that  this takes a very long time to run.
    clusterThresh = '3.1'; % cluster forming threshold. 
    nPerm       = '10000';% number of permutations
    % mask to run the permutation tests in and perform multiple comparisons in.
    % mask is intentionally of left hemi (lh) as the right hemi (rh) was
    % previously registered to the left, so it now has the lh indeces. This
    % happened inside searchlightDefinitionSurfaceWrapper.m
    mask          = fullfile(rootData,'freesurfer-subjects','fsaverage','label','lh.CORTEX.mgh');
    % run PALM
    groupLevelPermTests(rootData,'GLM2','xRun','surf',clusterThresh,nPerm,mask)
    
% # plot RDMs at the peaks: entorhinal for relational structure and 
%   LOC for visual identity  (both in right hemisphere)
%   First get vertex of peak activation, then plot average RDM across
%   participants in this vertex. 

% relational structure - Fig 2b and c, top
relationalStructMap = load_mgh(fullfile(rootData,'RDMs','GLM2','groupStats','perm','xRun','surf',...
                        ['relationalStructure_smth5_rh_nPerm10000_clstrTh3p1_dpv_tstat_uncp_c1.mgz']));
[~,indMax] = max(relationalStructMap);
roiStr = [num2str(indMax - 1),'rh']; % -1 is for conversion from Matlab to Freesurfer coords. 
plotDataRdm(rootData,subjects,'GLM2','xRun','relationalStructure',roiStr,fwhm)
plotGardnerAltman(rootData,subjects,'GLM2','xRun','surf','relationalStructure',roiStr,fwhm)

% visual identity - Fig 2b and c, bottom
visualIdentityMap = load_mgh(fullfile(rootData,'RDMs','GLM2','groupStats','perm','xRun','surf',...
                        ['visualIdentity_smth5_rh_nPerm10000_clstrTh3p1_dpv_tstat_uncp_c1.mgz']));
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

% run the RSA analysis on the cortical surface. 
% To perform statistical tests we will perform an ROI analysis on the peaks
% of the univatiate STRUCT prediction error effect from the same GLM (GLM3). 
for iSub=1:length(subjects)
    runRSA(rootData,subjects{iSub},'GLM3','xRun','surf',fwhm)
end

% # group analysis 
% first stack all single subjects contrasts into one file
stackSubjectsContrasts(rootData,'GLM3','xRun','surf');
% Run permutation tests in PALM.
% The permutation tests also correct for multiple comparisons. %
% However this takes a long time to run.
clusterThresh = 'None'; % cluster forming threshold. For this GLM we are 
                        % only looking at ROIs defined b y the peaks of the
                        % univarite prediction error effect. 
nPerm       = '10000';% number of permutations
% mask to run the permutation tests in and perform multiple comparisons in.
% mask is intentionally of left hemi (lh) as the right hemi (rh) was
% previously registered to the left, so it now has the lh indeces. This
% happened inside searchlightDefinitionSurfaceWrapper.m
mask          = fullfile(rootData,'freesurfer-subjects','fsaverage','label','lh.CORTEX.mgh');
% run PALM
groupLevelPermTests(rootData,'GLM3','xRun','surf',clusterThresh,nPerm,mask)

% plot subject-average data RDM in the vmPFC peak of univariate effect: Fig
% 3b and S5
plotDataRdm(rootData,subjects,'GLM3','xRun','peXstructure','51084lh',fwhm)
plotGardnerAltman(rootData,subjects,'GLM3','xRun','surf','peXstructure','51084lh',fwhm)

% # run a parallel pipeline in MNI space rather than subject space, in
% order to run RSA on volumetric searchlights (specifically in the ventral
% Striatim peak of the pe x structure interaction effect). 

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
stackSubjectsContrasts(rootData,'GLM3','xRun','lAccumbens20');
clusterThresh = 'None'; % searching only in a peak voxel (from the orthogonal univariate contrast).
nPerm       = '10000';% number of permutations
mask          = fullfile(rootData,'anatomical','mni152Masks','lAccumbens20.nii');
groupLevelPermTests(rootData,'GLM3','xRun','lAccumbens20',clusterThresh,nPerm,mask)

% Fig 3c: peak of ventral striatum effect is in MNI coords [-10,8,-12],
% corresponding to FSL voxel coords [50,67,30]
plotGardnerAltman(rootData,subjects,'GLM3','xRun','lAccumbens20','peXstructure',[50,67,30],fwhm)

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
stackSubjectsContrasts(rootData,'GLM2','corrWithinBlock','surf');
% Run permutation tests
% run on the surface permutation tests in PALM.
clusterThresh = '2.3'; % cluster forming threshold.
nPerm       = '10000';% number of permutations
% mask to run the permutation tests in and perform multiple comparisons in.
% correcMask is intentionally of left hemi (lh) as the right hemi (rh) was
% previously registered to the left, so it now has the lh indeces.
mask          = fullfile(rootData,'freesurfer-subjects','fsaverage','label','lh.CORTEX.mgh');
% run PALM
groupLevelPermTests(rootData,'GLM2','corrWithinBlock','surf',clusterThresh,nPerm,mask)

% plot RDMs at the peak of mOFC effect - Fig S6b. 
plotDataRdm(rootData,subjects,'GLM2','corrWithinBlock','latentVar_S-S','110371lh',fwhm)
plotDataRdm(rootData,subjects,'GLM2','corrWithinBlock','latentVar_O-O','110371lh',fwhm)
plotDataRdm(rootData,subjects,'GLM2','corrWithinBlock','latentVar_S-O','110371lh',fwhm)

% Fig 4a left
plotGardnerAltman(rootData,subjects,'GLM2','corrWithinBlock','surf','latentVar_all','110371lh',fwhm)
% fig 4b, left
plotGardnerAltman(rootData,subjects,'GLM2','corrWithinBlock','surf','latentVar_S-S','110371lh',fwhm)
% fig 4b, middle
plotGardnerAltman(rootData,subjects,'GLM2','corrWithinBlock','surf','latentVar_O-O','110371lh',fwhm)
% fig 4b, right
plotGardnerAltman(rootData,subjects,'GLM2','corrWithinBlock','surf','latentVar_S-O','110371lh',fwhm)

%% Subject-by-subject correaltion between peak vmPFC pe x structure & peak OFC latent variable effects - Fig S6c
plotVmpfcOfcEffectsCorr(rootData,subjects)

