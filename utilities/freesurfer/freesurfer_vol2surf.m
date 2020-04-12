function S=freesurfer_vol2surf(subj,images,varargin);
% Map multiple input volumes to left and right
% surface
% If no register.dat file is in the surface directory,
% a regsiter.dat is made
subj_dir=getenv('SUBJECTS_DIR');

hem={'lh','rh'};
surf={'.white'};
num_surf=5;
smoothing=1;
hemisphere=[1:2]; % Do both hemispheres

vararginoptions(varargin,{'smoothing'});


% Original subjects dir
surf_dir =[subj_dir filesep subj filesep 'surf'];
now_dir=pwd;
cd(surf_dir);


% If no register.dat exists, make one
if (~exist('register.dat','file'))
    freesurfer_mat2registerdat(subj,images{1});
end;

% Make sure the file does not have NaN's in it, replace with 0

for i=1:length(images)
    V=spm_vol(images{i});
    X=spm_read_vols(V);
    [a,name,ext,num]=spm_fileparts(images{i});
    if (any(isnan(X(:))))
        warning('file contains NaN, replacing with 0');
        X(isnan(X))=0;
        V.fname=[a filesep 'nn' name '.' ext]; 
        images{i}=V.fname;
        spm_write_vol(V,X);
    end;
    for h=hemisphere
        [stat,res]=system(['mri_vol2surf --mov ' images{i} ' --regheader ' subj ' --out '  hem{h} '.' name '.paint' ' --out_type paint --hemi ' hem{h}]);
        keyboard; 
    end;
end;
cd (now_dir);