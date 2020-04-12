function smoothVolumetricMapsWithinMask(rootData,inputDir,fnameFullPath,sl,fwhm)

% smooth an image (fname) within an appropriate mask - the one used to
% define sl. 

sigma  = num2str(fwhm/2.3548); % standard conversion 
mask   = fullfile(rootData,'anatomical','mni152Masks',[sl '.nii']);
fname_output = [fnameFullPath(1:end-4) '_smth' num2str(fwhm) '.nii.gz'];
cd(inputDir)

system(['fslmaths ' fnameFullPath ' -nan -mas ' mask ' -s ' sigma ' -mas ' mask ' result1']);
system(['fslmaths ' mask ' -s ' sigma ' -mas ' mask ' result2']);
system(['fslmaths result1 -div result2 ' fname_output]);
gunzip(fname_output);
delete('result*');
delete(fname_output);

