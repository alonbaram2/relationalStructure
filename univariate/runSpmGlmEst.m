function runSpmGlmEst(rootData,sub,glm)
% sub: subject name, e.g. sub-01
% glm: GLM name, e.g. GLM1

outputDir     = fullfile(rootData,'spmDirs',glm,sub);
evsDir        = fullfile(rootData,'evs',glm,sub); % folder where model spec files (i.e. evs/regressors) are, created in wrapper_univariate.m
preProcDir    = fullfile(rootData,'preprocessed',sub);
mcDir         = fullfile(rootData,'mc',sub); % folder wit motion correction files. in FSL this is produced during preproessing by MCFLIRT. 
anatomicalDir = fullfile(rootData,'anatomical',sub);
if ~exist(outputDir,'dir') % create a folder for 1st-level results
    mkdir(outputDir);
end

% aquisition parameters
TR      = 1.235;
nSlices = 72;

% mask for running GLMs and getting the right dimensions. I moved all the 
% functional data to 2x2x2 version of the structural space (created with 
%FLIRT, usin FSL's registration tools), because the
% structural scan was taken in the middle of the session (after 4 EPI
% blocks, before another 4 EPI blocks). This step is not necessary - data
% can be in the original functional space (as long as different EPI blocks
% are aligned) and then the mask here should be a liberal mask in that
% space. 
mask = fullfile(anatomicalDir,'structural_brain_mask_2mm.nii');

blockNames = {'run1_cond1','run1_cond2','run1_cond3','run1_cond4',...
    'run2_cond1','run2_cond2','run2_cond3','run2_cond4'}; % in the SPM nomenclature what I call "blocks" 
% will be called "runs" - i.e. different EPI sequences. 

%% set up SPM batch
% J is what is ususally called "fmri_spec" in SPM code
J.dir = {outputDir};
J.timing.units = 'secs';
J.timing.RT = TR;
J.timing.fmri_t = nSlices; % number of slices
J.timing.fmri_t0 = 1;        % reference slice
for iBlock=1:numel(blockNames) % for each block/
    [filePathListThisRun] = spm_select('FPList',fullfile(preProcDir),[blockNames{iBlock},'^*.*nii$']);
    filePathListThisRun = mat2cell(filePathListThisRun, ones(size(filePathListThisRun,1),1));%convert to cell arrays
    % scans
    J.sess(iBlock).scans=filePathListThisRun;
    % set to empty
    J.sess(iBlock).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
    % multi cond files
    J.sess(iBlock).multi = {fullfile(evsDir,[blockNames{iBlock},'.mat'])}; 
    % set to empty
    J.sess(iBlock).regress = struct('name', {}, 'val', {});
        
    % transform motion parameter files outputed by FSL preprocessing to a
    % format that SPM can understand. Note that the order of the output of FSL preprocessingis different 
    % to the usual order in SPM - in FSL it's angles (in rad) first (roll,
    % yaw, pitch), then translation (x,y,z). 
    R = load(fullfile(mcDir,[blockNames{iBlock},'.par'])); % load motion correction parameters.
    names = {'roll','yaw','pitch','x','y','z'}; % names of columns in motion parameters files of FSL. 
    save(fullfile(mcDir,[blockNames{iBlock},'.mat']),'names','R');    
    J.sess(iBlock).multi_reg = {fullfile(mcDir,[blockNames{iBlock},'.mat'])};      
    % high pass filter
    J.sess(iBlock).hpf = inf; % I already high-passed filtered during preprocessing using SPM, so no need to filter here. I will filter the design matrix manually later. 
end
J.fact = struct('name', {}, 'levels', {});
J.bases.hrf.derivs = [0 0];
J.volt = 1;
J.global = 'None';
J.mthresh = 0.000001;  % we will use an explicit mask of the brain, so setting this to very low. 
J.mask = {[mask,',1']};
J.cvi = 'fast';   % autocorrelation correction
matlabbatch{1}.spm.stats.fmri_spec = J;
% model estimation
estName=fullfile(outputDir,'SPM.mat'); % SPM.mat for estimation
matlabbatch{2}.spm.stats.fmri_est.spmmat = {estName};
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
% run model specification and estation
spm_jobman('run',matlabbatch);