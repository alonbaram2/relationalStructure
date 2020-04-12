function runRSA_withinCluster(rootData,sub,glm,distType,sl,fwhm,iBlock)
% This is adapted from the RSA toolbox, see
% https://github.com/rsagroup/rsatoolbox

% distype: xRun (cross-run correlation) or corr (regular correaltion
% distance)
% sl: searchlight (surf or striatum)


mkdir(fullfile(rootData,'RDMs',glm,sub))
cd(fullfile(rootData,'RDMs',glm,sub))

% load seachlight definitions (surface)
if strcmp(sl,'surf')
    L = load(fullfile(rootData,'searchlights',sub,['SL_' sl '_v100.mat']));
    % load SPM structures of GLMs in subject space
    load(fullfile(rootData,'spmDirs',glm,sub,'SPM.mat'));
else % volumetric searchlight
    L = load(fullfile(rootData,'searchlights',sub,['SL_vol_' sl '_v100.mat']));
    % load SPM structures of GLMs in MNI space
    load(fullfile(rootData,'spmDirs_MNI',glm,sub,'SPM.mat'));
end

switch distType
    case 'xRun'
 
        
    case 'corrWithinBlock'
        % within-block latent variable representation: Fig 4.
        % This also uses GLM2, but it uses the betas for both stimulus and
        % outcome presentation time for all three stimuli. The distances
        % are also calculated differently - first the within-block
        % distances (simple correlation distance) are calculated, and then
        % the RDMs of the different blocks are averaged.
        
%         for iBlock =1:8 % this loop can easily be parallelised
            % Get all the regressors that take a part in the RSA analysis:
            % logical vector the length of all regressors
            condition = false(1,size(SPM.xX.xKXs.X,2));
            regInds4RSA = [1,2,3,4,5,6]; % 1,2,3: [A], [B], [C] at stim presentation. 4,5,6,: [A], [B], [C] at outcome presentation
            condition(SPM.Sess(iBlock).col(regInds4RSA)) = true;
            % Get number of elements in within-block RDM.
            nElem = sum(condition);
            
            % file name for RDM
            outFiles={};
            for k=1:nElem*(nElem-1)/2 % k is 4th dim of RDM NIFTI.
                outFiles{end+1}=sprintf('dist_%s_%s_bl%d.nii,%d',distType,sl,iBlock,k);
            end;
            
            analysisName =[glm,'_',distType,'_',sl,'_bl',int2str(iBlock)];
            % run the analysis
            rsa.runSearchlight(L,SPM.xY.VY,outFiles,analysisName,@calcRDMs,'optionalParams',{SPM,condition});

            
%         end
end



