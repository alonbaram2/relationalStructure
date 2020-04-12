function caret_vol2surf(spec,images,varargin);
% function caret_vol2surf(spec,images,varargin);
% Map multiple input volumes to a spec file (FIDUCIAL surface)
% Adds the file to the spec file 
% Old version of the file is overwritten, unless you specify mode to be 'append' 
% INPUT: 
%   spec: Path and filename of spec-file
%   images: cell array of images to be mapped (these need to be prepared
%           with caret_volprep
%   
% VARARGIN 
%   'algorithm',mapping algorithm (default: METRIC_GAUSSIAN)
%           other choices:  METRIC_AVERAGE_VOXE
%                           METRIC_ENCLOSING_VOXEL
%                           PAINT_ENCLOSING_VOXEL
%   'gaussparam',[neighborhood ,...
%                 sigma_normal,...
%                 sigma_tangent,... 
%                 normal_below_cutoff,...
%                 normal_above_cutoff,... 
%                 tangent_cutoff]
%    'cubesize',2  : parameter for METRIC_AVERAGE_VOXEL   
% ---------------------------
% v1.0 Joern Diedrichsen (j.diedrichsen@ucl.ac.uk)
%
% v.1.1: Tobias Wiestler (31.05.2011): Force to write binary metric files 
% v.1.2: Joern added informative error message when caret_command fails
% 
% POSSIBLE EXTENSION: USE BOTH PIAL AND WHITE SURFACE TO DO MAPPING

[sdir,specname]=fileparts(spec); 
SPEC=caret_load(spec); 

% Settings 
algorithm='METRIC_AVERAGE_VOXEL';
algorithm='METRIC_GAUSSIAN'; 
gaussparam=[4 2 2 1 3 3]; 
cubesize=2; 
mode='replace'; 
outname='func.metric'; 

vararginoptions(varargin,{'smoothing','outname','mode','algorithm'});
now_dir=pwd; 
cd(sdir);
filelist=[];

% Make sure the file does not have NaN's in it, replace with 0
for i=1:length(images)
    V=spm_vol(images{i});
    filelist=[filelist ' ' images{i}];
end;

switch (mode) 
    case 'replace'
        inname='""'; 
    case 'append' 
        inname=outname;
    otherwise
        error('Invalid mode: replace / append'); 
end; 

comm=['caret_command -volume-map-to-surface ' ...
        SPEC.FIDUCIALcoord_file ' '  SPEC.CLOSEDtopo_file ' ' ...
        inname ' ' outname ' ' algorithm ' ' filelist ' -WRITE-FILE-FORMAT-METRIC BINARY']; 
    
switch(algorithm) 
    case 'METRIC_GAUSSIAN' 
        comm=[comm ' -g ' num2str(gaussparam,'%2.2f ')]; 
        
    case 'METRIC_AVERAGE_VOXEL'  
        comm=[comm ' -av ' num2str(cubesize)]; 
        
end; 
[stat,res]=system(comm);
if (stat~=0)
    error(['Error during caret-command: ' res]);
end; 
if (~isempty(SPEC.data{end}))
    SPEC.data{end+1}=''; 
end; 
SPEC.data{end}=['metric_file ' outname]; 
caret_save(spec,SPEC); 
cd (now_dir);