function runRSA(rootData,sub,glm,distType,sl,fwhm)

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
        analysisName =[glm,'_',distType,'_',sl];
        
        % Get all the regressors that take a part in the RSA analysis:
        % logical vector the length of all regressors
        condition = false(1,size(SPM.xX.xKXs.X,2));
        for iBlock =1:8
            % Indeces  of relevant regressors in a specific block, see getEvs.m
            switch glm
                case 'GLM2'  % Visual identity and Relational structure between stimuli: Fig 2.
                    % [A], [B], [C] at stim presentation. note that we need
                    % Note that we need the [C] trials elements for the visual
                    % identity analysis. We will ignore them for the
                    % relational structure analysis.
                    regInds4RSA = [1,2,3];
                case 'GLM3'
                    regInds4RSA = 2; % crPe parmaetric modulation on [AB] trials at outcome time
            end
            condition(SPM.Sess(iBlock).col(regInds4RSA)) = true;
        end
        % number of elements in the exapnded RDM, before averaging across runs.
        % Note that here we find both within and across-run distances, but we will
        % only use the across runs distances (this is an unefficient way of
        % computing this but it's easier to undrstand what's going on).
        nElemBeforeCollapse = sum(condition);
        
        % Length of RDM after collapsing the two across-run distances for each
        % couple of conditions. This is the length of the final RDM, as will
        % appear in the plots (though for relational structure analysis we will
        % ifnor C trials). It is also the number of different experimental
        % conditions.
        nElemFinal = nElemBeforeCollapse/2;
        
        % file name for RDM
        outFiles={};
        for k=1:nElemBeforeCollapse*(nElemBeforeCollapse-1)/2 % k is 4th dim of RDM NIFTI.
            outFiles{end+1}=sprintf('dist_%s_%s.nii,%d',distType,sl,k);
        end;
        
        % run the analysis
        rsa.runSearchlight(L,SPM.xY.VY,outFiles,analysisName,@calcRDMs,'optionalParams',{SPM,condition});
        
        % # Average the two cross-run distances between every couple of
        % conditions. Diagonal of RDM is saved in the last nElemFinal elemensts.
        % This will write over the expanded (nElemBeforeCollapse x nElemBeforeCollapse)
        % rdm created in the last line, as we don't need it anymore.
        collapseRdmXRun(rootData,sub,glm,sl,nElemFinal)
        
        % calculate a summary statistic per voxel per subject, project to the surface 
        % if applicable) and smooth (within the relevant searchlight/on the
        % surface). 
        withinSubjRsaStats(rootData,sub,glm,distType,sl,fwhm)   
        

        
    case 'corrWithinBlock'
        % within-block latent variable representation: Fig 4.
        % This also uses GLM2, but it uses the betas for both stimulus and
        % outcome presentation time for all three stimuli. The distances
        % are also calculated differently - first the within-block
        % distances (simple correlation distance) are calculated, and then
        % the RDMs of the different blocks are averaged.
        
        for iBlock =1:8 % this loop can easily be parallelised
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
        end
        
        % average RDMs from the different blocks.
        system(['fslmaths dist_' distType '_' sl '_bl1 -add dist_' distType '_' sl '_bl2 -add ' ...
            'dist_' distType '_' sl '_bl3 -add dist_' distType '_' sl '_bl4 -add ' ...
            'dist_' distType '_' sl '_bl5 -add dist_' distType '_' sl '_bl6 -add ' ...
            'dist_' distType '_' sl '_bl7 -add dist_' distType '_' sl '_bl8 ' ...
            '-div 8 dist_' distType '_' sl '_mean']);
        gunzip(['dist_' distType '_' sl '_mean.nii.gz'])
        delete(['dist_' distType '_' sl '_mean.nii.gz'])
        
        % calculate a summary statistic per voxel per subject, project to the surface 
        % if applicable) and smooth (within the relevant searchlight/on the
        % surface). 
        withinSubjRsaStats(rootData,sub,glm,distType,sl,fwhm)        

end



