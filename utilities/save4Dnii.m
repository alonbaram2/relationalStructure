
function save4Dnii(nii,V,filename)

% nii: 4D or 3D matrix to save.
% V: either an SPM structure or just the output of spm_vols - i.e a
%    structure with NIFTI header information.

if ndims(nii) == 3
    T = 1; % no 4th dimension
elseif ndims(nii) == 4
    T = size(nii,4);
end
if isempty(V)
    dim = [size(nii,1),size(nii,2),size(nii,3)];
    mat = [-2,0,0,92 ; 0,2,0,-128 ; 0,0,2,-74 ; 0,0,0,1]; % this is taken from an MNI 152 2mm^3. 
    warning('orientation of saved nifti file might be incorrect - consider using fslcpgeom with flag -d to copy header information from a NIFTI file with the same desired orientation.')
else
    if isfield(V,'xY') % V is an SPM structure
        dim = V.xY.VY(1).dim;
        mat = V.xY.VY(1).mat;
    else
        dim = V(1).dim;
        mat = V(1).mat;
    end
end
for i=1:T
    Z = squeeze(nii(:,:,:,i));
    Vo      = struct(...
        'fname',    fullfile(filename),...
        'dim',      dim,...
        'dt',       [spm_type('float32') spm_platform('bigend')],...
        'mat',      mat,...
        'n',        [i 1],...
        'descrip',  '');
    spm_write_vol(Vo,Z);
end