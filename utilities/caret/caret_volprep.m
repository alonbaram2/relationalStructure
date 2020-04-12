function images=caret_volprep(images,varargin);
% images=caret_volperp(images,varargin); 
% Takes a list of file names and prepares them 
% for surface mapping under caret: reslices them into LPI and removes possible NaN's. 
% if no outnames are given, the new volumes are saved as ['r' name '.nii'] 
% INPUT: 
%   images: cell array of images to be prepared 
% VARARGIN: 
%   'vol_padding',val:  replace NaN with val (default: 0)
%   'outnames',{}:      cell array of output files names 
% ---------------------------
% v1.0 Joern Diedrichsen (j.diedrichsen@ucl.ac.uk
% 
vol_padding=0;      % Volume padding 
outnames={}; 
vararginoptions(varargin,{'vol_padding','outnames'});

for i=1:length(images)
    [path,name,ext,num]=spm_fileparts(images{i});
    if (isempty(outnames))
        outname=fullfile(path,['r' name  '.nii' num]);
    else 
        outname=outnames{i};
    end; 
    images{i}=spmj_reslice_LPI(images{i},'name',outname); 
    
    % Check for NaN
    V=spm_vol(images{i}); 
    X=spm_read_vols(V); 
    if (any(isnan(X(:))))
        X(isnan(X))=vol_padding; 
        spm_write_vol(V,X); 
    end; 
end; 
