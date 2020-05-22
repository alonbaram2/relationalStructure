function projectTMapsToSurface(rootData,glm)
    
tmp = dir(fullfile(rootData,'spmDirs',glm,'2ndLevel'));
conName = tmp(3).name; % this assumes there's only one contrast to project - if there are more it will need adjusting.
conDir = fullfile(rootData,'spmDirs',glm,'2ndLevel',conName);
cd(conDir);
fname = dir('spmT*.nii');
fname = fname.name(1:end-4);
 
% project to surface, smooth with 1mm on surface for visualisation
system(['DYLD_LIBRARY_PATH=' getenv('FREESURFER_HOME') '/lib/gcc/lib mri_vol2surf --trgsubject  fsaverage --mov ' fname '.nii --o ' fname '_lh.mgh  --reg ' getenv('FREESURFER_HOME') '/average/mni152.register.dat --cortex --fwhm 1 --hemi lh'])
system(['DYLD_LIBRARY_PATH=' getenv('FREESURFER_HOME') '/lib/gcc/lib mri_vol2surf --trgsubject  fsaverage --mov ' fname '.nii --o ' fname '_rh.mgh  --reg ' getenv('FREESURFER_HOME') '/average/mni152.register.dat --cortex --fwhm 2 --hemi rh'])

% flip rh hemi to lh indeces - this is for compatibility with the rest of
% the paper where the rh is flipped (and should be presented on the lh
% surface)
system(['DYLD_LIBRARY_PATH=' getenv('FREESURFER_HOME') '/lib/gcc/lib mris_apply_reg --src ' fname '_rh.mgh --trg ' fname '_rh.mgh --streg ' getenv('SUBJECTS_DIR') 'fsaverage/surf/lh.sphere.left_right ' getenv('SUBJECTS_DIR') '/fsaverage/surf/rh.sphere.left_right']);