function searchlightDefinitionSurfaceWrapper(rootData, sub)
% Define searchlights on subject's surface.
% This code assumes you have already ran the Freesurfer surface reconstruciton
% using the recon-all command: https://surfer.nmr.mgh.harvard.edu/fswiki/recon-all
% This code was written by Joern Diedrichsen (j.diedrichsen@ucl.ac.uk) and 
% was taken from the RSA toolbox github (https://github.com/rsagroup/rsatoolbox). 

% This code uses the Caret software (now depricated for interface between
% Freesurfer and Matlab. 
% These type of analyses are nowadays done using WorkBench (see https://osf.io/k89fh/wiki/Surface/)
% though I did not use it in the analysis for this paper.

% I've added comments to explain what happens at each stage. 

%%%%%%%%

% get the Freesurfer SUBJECTS_DIR. This assumes MATLAB inherits the environment
% variables from the shell, and that you've set up the FREESURFER variables.
subject_dir = getenv('SUBJECTS_DIR');
anatomicalDir = fullfile(rootData,'anatomical',sub);

if exist(fullfile(anatomicalDir,'structural_brain_mask_2mm.nii.gz'),'file')
    gunzip(fullfile(anatomicalDir,'structural_brain_mask_2mm.nii.gz'));
end
% mask to find searchlights in. Only the voxels corresponding to the
% cortical surface within this mask will be used. 
Vmask = rsa.readMask(fullfile(anatomicalDir,'structural_brain_mask_2mm.nii')); % mask
% File name of searchlight matlab structure
fname = 'SL_surf_v100.mat'; % 100 voxels in each searchlight
% make dir to tore searchlight file
slDir = fullfile(rootData,'searchlights',sub);
if ~exist(slDir,'dir') % create a folder for 1st-level results
    mkdir(slDir);
end
% make dir to store Caret files - surface files resampled 
caretDir = fullfile(rootData,'caretDirs',sub);
if ~exist(caretDir,'dir') % create a folder for 1st-level results
    mkdir(caretDir);
end
if ~exist(fullfile(subject_dir,['x' sub]) ,'dir')
    
    % register the right hemisphere to the left (xhemi) and register
    % both to fsaverage_sym (atlas), using Freesurfer commands.
    % First inverts the right hemi using 'xhemireg'.
    % Then creates registration files using 'surfreg':
    % Left hemi: lh.fsaverage_sym.sphere.reg in subject/surf
    % Right Hemi: lh.fsaverage_sym.sphere.reg in subject/xhemi/surf
    % The reason for doing this is convenience: it will mean that the same
    % indexing of verteces will apply for both left and right hemisphere,
    % such that  vertex i in both hemispheres will be in the same
    % anatomical location.
    freesurfer_registerXhem({sub},subject_dir,'hemisphere',[1 2]);
    
    % resample the registered subject's surfaces into atlas space (ico7). Creates new
    % remapped surfaces (.white, .pial, .inflated) and curves (.cure, .sulc, .area)
    % which will be in $SUBJECTS_DIR/['x' subj] dir (both hemisphere).
    % Effectively this means changeing the number and order of verteces of the input
    % surfaces to match ico7, but keep the information about the
    % the RAS coordinates of the nearest neighbour vertex
    % in the (moved) subject surface to each ico vertex.
    % the topology (surface faces) is taken from ico.
    freesurfer_mapicosahedron_xhem(sub,subject_dir,'smoothing',1,'hemisphere',[1:2]);
    
    % Move the surfaces vertex coords to correspond to World
    % rather than RAS_tkreg coords (256x256x256,
    % like in $SUBJECTS_DIR/$sub/orig.mgz).
    % save in Caret format.
    caret_importfreesurfer(['x',sub],subject_dir,caretDir);
end
LcaretDir = fullfile(caretDir,sprintf(['x',sub]),'LeftHem');
RcaretDir = fullfile(caretDir,sprintf(['x',sub]),'RightHem');
white     = {fullfile(LcaretDir,'lh.WHITE.coord'),fullfile(RcaretDir,'rh.WHITE.coord')};
pial      = {fullfile(LcaretDir,'lh.PIAL.coord'),fullfile(RcaretDir,'rh.PIAL.coord')};
topo      = {fullfile(LcaretDir,'lh.CLOSED.topo'),fullfile(RcaretDir,'rh.CLOSED.topo')};

% S is a structure with the subjects surfaces in world
% (mm) coords. I.e. if you save S as a surface and load it in
% Freeview, the "RAS" coordinates of each vertex will correspond
% exactly to the "scanner anatomical" coords in FSLeyes.
S         = rsa.readSurf(white,pial,topo);
L = rsa.defineSearchlight_surface(S,Vmask,'sphere',[10, 100]); % searchloght of 100 voxels, max radius 10mm
save(fullfile(slDir,fname),'-struct','L'); 
