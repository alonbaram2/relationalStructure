function plotGardnerAltman(rootData,subjects,glm,distType,sl,analysisName,roi,fwhm)

% plot using "dabest" package for plotting diffeerences between two groups:
% https://www.mathworks.com/matlabcentral/fileexchange/65260-dabest-matlab

% if sl==surf, roi is a sreing of the vertex number followed by hemi: e.g. 137811rh. 
% otherwise, roi is a vector of the FSLeyes coordinates of the requested voxel. 
% e.g [30,68,10]. Note that these coordinates are 0-based, and are different.
% from the MNI coordinates. 

similar = nan(length(subjects),1);
dissimilar = nan(length(subjects),1);

for iSub=1:length(subjects)
    sub=subjects{iSub};        
    inDir = fullfile(rootData,'RDMs',glm,sub,'rsaContrastsGroups',distType,sl);
    if strcmp(sl,'surf')
    dataFile_similar = fullfile(inDir,[analysisName,'_sim_smth' int2str(fwhm) '_' roi(end-1:end) '.mgh']);
    dataFile_dissimilar = fullfile(inDir,[analysisName,'_dsim_smth' int2str(fwhm) '_' roi(end-1:end) '.mgh']);
    tmpSim = load_mgh(dataFile_similar);
    similar(iSub) = tmpSim(str2num(roi(1:end-2)) +1); % add 1 for conversion from 0-based freesurfer to 1-based Matlab
    tmpDissim = load_mgh(dataFile_dissimilar);
    dissimilar(iSub) = tmpDissim(str2num(roi(1:end-2)) +1); % add 1 for conversion from 0-based freesurfer to 1-based Matlab

    else
        dataFile_similar = fullfile(inDir,[analysisName,'_sim_smth' int2str(fwhm) '.nii']);
        dataFile_dissimilar = fullfile(inDir,[analysisName,'_dsim_smth' int2str(fwhm) '.nii']);
        tmpSim = read_avw(dataFile_similar);
        similar(iSub) = tmpSim(roi(1)+1,roi(2)+1,roi(3)+1); % add 1 for conversion from 0-based freesurfer to 1-based Matlab
        tmpDissim = read_avw(dataFile_dissimilar);
        dissimilar(iSub) = tmpDissim(roi(1)+1,roi(2)+1,roi(3)+1); % add 1 for conversion from 0-based freesurfer to 1-based Matlab
    end
end

% plot

identifiers = cell(2*length(similar),1);
identifiers(1:length(similar)) = {'Same Struct'};
identifiers(1+length(similar):end) = {'Diff Struct'};
data = vertcat(similar,dissimilar);
dabest(identifiers,data,'Paired');