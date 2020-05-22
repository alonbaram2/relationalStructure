function stackSubjectsContrasts(rootData,subjects,glm,distType,sl)

% get analyses names from example subject
exampleInputDir = fullfile(rootData,'RDMs',glm,'sub-01','rsaContrasts',distType,sl);
cd(exampleInputDir)
tmpFname = dir('*smth*'); % find all files that were already smoothed - they are the ones to use for concatination
fname = cell(size(tmpFname));
for iAnalysis = 1:length(fname)
    fname{iAnalysis} = tmpFname(iAnalysis).name(1:end-4); % get rid of .nii suffix
end

outputDir = fullfile(rootData,'RDMs',glm,'groupStats','stackedInputs',distType,sl);
mkdir(outputDir);

for iAn = 1:length(fname)
    if strcmp(sl,'surf')
            str = ['mri_concat '];
            for iSub = 1:length(subjects)
                str = [str fullfile(rootData,'RDMs',glm,subjects{iSub},'rsaContrasts',distType,sl,[fname{iAn} '.mgh']) ' '];
            end
            str = [str '--o ' fullfile(outputDir,[fname{iAn} '.mgh'])];
            system (str);           
    else            
            str = ['fslmerge -t ' fullfile(outputDir,fname{iAn}) ' '];
            for iSub = 1:length(subjects)
                str = [str fullfile(rootData,'RDMs',glm,subjects{iSub},'rsaContrasts',distType,sl,[fname{iAn} '.nii']) ' '];
            end
            system (str);
            gunzip(fullfile(outputDir,[fname{iAn} '.nii.gz']));
            delete(fullfile(outputDir,[fname{iAn} '.nii.gz']))
    end    
end