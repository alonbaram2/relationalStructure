% Example function to demonstrate surface alignment for surface-based voxel
% selection

% set directory
d='/Volumes/organized/116_surfing_toolbox/fingerdata-0.2/';

% set up variables
C=struct();
C.action='all';
C.subjid='s88';
C.origvolfn=[d 'glm/anat_al+orig'];
C.fssurfdir=[d 'fs/s88/surf/'];
C.aligndir=[d 'ref/'];

% run alignment
[cmd,S]=surfing_suma_align(C);

% make SUMA spec files
surfing_suma_makespec(S);

% hardlink to glm directory
srcdir=C.aligndir;
trgdir=[d 'glm/'];
surfing_suma_surfacefiles(srcdir,trgdir,'ln');