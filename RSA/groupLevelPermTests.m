function groupLevelPermTests(rootData,glm,distType,sl,clusterThresh,nPerm,maskName)

% make sure palm is installed an in path.
palmDir = '/home/fs0/abaram/scratch/MATLAB/palm-alpha111/';
addpath(genpath(palmDir));

groupStatsStr = 'groupStats_all'; % groupStats_all: all 28 subjects; groupStats: just the 25 that consented to share data

% mask to run the permutation tests in and perform multiple comparisons in.
% mask is intentionally of left hemi (lh) as the right hemi (rh) was
% previously registered to the left, so it now has the lh indeces. This
% happened inside searchlightDefinitionSurfaceWrapper.m
maskFile = fullfile(rootData,'masksAndRois','fsaverage',['lh.' maskName '.mgh']);

inDir   = fullfile(rootData,'RDMs',glm,groupStatsStr,'stackedInputs',distType,sl);
cd(inDir)
fname = dir('*smth*'); % find all files that were already smoothed - they are the ones to use for concatination

outDir = fullfile(rootData,'RDMs',glm,groupStatsStr,'perm',distType,sl,maskName);
mkdir(outDir);

if strcmp(sl,'surf')
    surface = 'pial';
    surfDir  = fullfile(rootData,'freesurfer-subjects','fsaverage','surf');
    surfFile = fullfile([surfDir '/lh.' surface]); % this is intentionally always lh because rh data was flipped to be with lh indeces.
    
    if ~exist(maskFile,'file')
        maskLabelFName = [maskFile(1:end-4) '.label'];
        % change format of mask from .label to .mgh as that's what PALM
        % wants
        system (['mri_label2label --s fsaverage --srclabel ' maskLabelFName ...
            ' --trglabel ' maskLabelFName ' --regmethod surface --hemi lh --outmask ' maskFile]);
    end
    for iAnalysis = 1:length(fname) % better to parallelise this
        outFile = fullfile(outDir,[fname(iAnalysis).name(1:end-4) '_nPerm' nPerm '_clstrTh' clusterThresh]);        
        outFile = strrep(outFile,'.','p'); % switch decimal point to 'p', as PALM doesn't like points in filenames.        
        str = ['palm -i ' fname(iAnalysis).name ' -s ' surfFile ...
            ' -n ' nPerm ' -o ' outFile ' -ise -save1-p -m ' maskFile];
        if ~strcmp(clusterThresh,'None')
            str = [str ' -C ' clusterThresh ' -Cstat mass' ];
        end
        eval(str)
    end    
else % volumetric
    for iAnalysis = 1:length(fname) % better to parallelise this
        outFile = fullfile(outDir,[fname(iAnalysis).name(1:end-4) '_nPerm' nPerm '_clstrTh' clusterThresh]);
        outFile = strrep(outFile,'.','p'); % switch decimal point to 'p', as PALM doesn't like points in filenames.
        str = ['palm -i ' fname(iAnalysis).name ' -n ' nPerm ' -o ' outFile ' -ise -save1-p -m ' maskFile];
        if ~strcmp(clusterThresh,'None')
            str = [str ' -C ' clusterThresh ' -Cstat mass' ];
        end        
        eval(str)
    end
end


