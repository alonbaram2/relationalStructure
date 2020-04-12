function runSpmGlmEst_MNI(rootData,sub,glm)
% This function is identical to runSpmGlmEst, only it assumes the
% preprocessed data is in MNI space, saves the outputs in appropriate
% directories (with _MNI suffix), and uses the MNI mask.

% sub: subject name, e.g. sub-01
% glm: GLM name, e.g. GLM1

outputDir     = fullfile(rootData,'spmDirs_MNI',glm,sub);
evsDir        = fullfile(rootData,'evs',glm,sub); % folder where model spec files (i.e. evs/regressors) are, created in wrapper_univariate.m
preProcDir    = fullfile(rootData,'preprocessed_MNI',sub);
mcDir         = fullfile(rootData,'mc',sub); % folder wit motion correction files. in FSL this is produced during preproessing by MCFLIRT. 
anatomicalDir = fullfile(rootData,'anatomical',sub);
if ~exist(outputDir,'dir') % create a folder for 1st-level results
    mkdir(outputDir);
end

% aquisition parameters
TR      = 1.235;
nSlices = 72;

mask = fullfile(rootData,'anatomical','mni152Masks','brain.nii'); % in MNI space! taken from FSL standard

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
    % format that SPM can understand 
    R = load(fullfile(mcDir,[blockNames{iBlock},'.par'])); % load motion correction parameters.
    names = {'x','y','z','roll','yaw','pitch'}; % names of columns in motion parameters files
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