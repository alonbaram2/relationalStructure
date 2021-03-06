function withinSubjRsaStats(rootData,sub,glm,distType,sl,fwhm)

saveContrastsFlag        = true; % this is to get the whole brain group level maps
saveRdmElementsMapsFlag  = true; % this will be needed get the data RDM at specific ROis
% we need to save the elements as
% different files as we need to separately
% smooth each element on the surface.

rdmDir    = fullfile(rootData,'RDMs',glm,sub);

switch distType
    case 'xRun'
        
        % load data
        dataFile = fullfile(rdmDir,['dist_xRun_' sl '_collapsed.nii']);
        nii      = read_avw(dataFile);
        V        = spm_vol(dataFile); % headre info
        
        % Names of the analyses to run. These are the different model RDMs
        % that will be uset to calculate summary statistics from the data
        % RDM.
        switch glm
            case 'GLM2'
                % This are all the analyses for figure 2 and S2.
                analyses.names = {'visualIdentity','relationalStructure','relationalStructure_control1','relationalStructure_control2'};
            case 'GLM3'
                analyses.names = {'peXstructure'};
        end
    case 'corrWithinBlock'
        % GLM2 within block analyses
        
        % load data - mean RDM across block
        dataFile = fullfile(rdmDir,['dist_corrWithinBlock_' sl '_mean.nii']);
        nii      = read_avw(dataFile);
        V        = spm_vol(dataFile); % header info
        
        analyses.names = {'latentVar_all','latentVar_S-S','latentVar_O-O','latentVar_S-O'};
end
% get contrasts - this is to get the whole brain maps.
if saveContrastsFlag
    for iAn = 1: length(analyses.names)
        % get model RDM
        analyses.RDMs{iAn} = getModelRdm(analyses.names{iAn},false);
        
        outputDir = fullfile(rootData,'RDMs',glm,sub,'rsaContrasts',distType,sl);
        mkdir(outputDir);
        cd (outputDir)
                       
        % calculate summary statitic
        dissimilarMean = squeeze(nanmean(nii(:,:,:,analyses.RDMs{iAn}==1),4));
        similarMean = squeeze(nanmean(nii(:,:,:,analyses.RDMs{iAn}==0),4));
        analyses.statistic{iAn} = dissimilarMean - similarMean;
        
        % save contrast nifti. 
        fname_con  = [analyses.names{iAn} '.nii'];        
        save4Dnii(analyses.statistic{iAn},V,fullfile(outputDir,fname_con))                              
                
        % project to the surface if applicable, and smooth on the surface.
        % Otherwise smooth within the volumetric mask used to create the
        % searchlight.
        if strcmp(sl,'surf') % surface analyses
            %  project to the surface, and smooth on the surface
            projectToSurfAndSmooth(rootData,outputDir,fname_con,sub,fwhm);
        else % volumetric analyses - smooth within volumetric mask
            smoothVolumetricMapsWithinMask(rootData,outputDir,fullfile(outputDir,fname_con),sl,fwhm);
        end
        
        % do the same for the simialr and dissimilar groups of elements.
        % we'll need those to produce the gGardner-Altman plots.
        outputDir = fullfile(rootData,'RDMs',glm,sub,'rsaContrastsGroups',distType,sl);
        mkdir(outputDir);
        cd (outputDir)        
        fname_sim  = [analyses.names{iAn} '_sim.nii'];
        fname_dsim = [analyses.names{iAn} '_dsim.nii'];
        save4Dnii(similarMean,V,fullfile(outputDir,fname_sim))
        save4Dnii(dissimilarMean,V,fullfile(outputDir,fname_dsim))
        if strcmp(sl,'surf') % surface analyses
            projectToSurfAndSmooth(rootData,outputDir,fname_sim,sub,fwhm);
            projectToSurfAndSmooth(rootData,outputDir,fname_dsim,sub,fwhm);
        else % volumetric analyses - smooth within volumetric mask
            smoothVolumetricMapsWithinMask(rootData,outputDir,fullfile(outputDir,fname_sim),sl,fwhm);
            smoothVolumetricMapsWithinMask(rootData,outputDir,fullfile(outputDir,fname_dsim),sl,fwhm);            
        end                       
    end
end
% save individual RDM elements as separate files
% this will be needed get the data RDM at specific ROis
if saveRdmElementsMapsFlag
    outputDir = fullfile(rootData,'RDMs',glm,sub,'rdmElements',distType,sl);
    mkdir(outputDir);
    cd (outputDir)
    
    for iElem=1:size(nii,4)
        elemMap = squeeze(nii(:,:,:,iElem));
        fname_elem = ['elem' num2str(iElem) '.nii'];
        save4Dnii(elemMap,V,fullfile(outputDir,fname_elem))
        if strcmp(sl,'surf') % surface analyses
            %  project to the surface, and smooth on the surface
            projectToSurfAndSmooth(rootData,outputDir,fname_elem,sub,fwhm);
        else % volumetric analyses - smooth within volumetric mask
            smoothVolumetricMapsWithinMask(rootData,outputDir,fullfile(outputDir,fname_elem),sl,fwhm);
        end
    end
end
