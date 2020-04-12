function projectToSurfAndSmooth(rootData,outputDir,fname,sub,fwhm)
% note we have already changed directory to outputDir before calling this
% function. 

hemName = {'LeftHem','RightHem'};
hem     = {'lh','rh'};
caretDir = fullfile(rootData,'caretDirs',sub);
[~, faces] = read_surf(fullfile(rootData,'freesurfer-subjects','fsaverage','surf','lh.sphere.reg'));

for h=1:2
    % the surfaces in caretSDir are already in world (not
    % RAS) coordinates
    caretSDir = fullfile(caretDir,sprintf(['x',sub]),hemName{h});
    white     = caret_load(fullfile(caretSDir,[hem{h} '.WHITE.coord']));
    pial      = caret_load(fullfile(caretSDir,[hem{h} '.PIAL.coord']));
    images    = fullfile(outputDir,fname);
    outfile   = [fname(1:end-4) '.metric'];
    % map data onto surface, according to world coords.
    % notice vertex indeces correspond (for both hemi) to the left hemi of
    % fsaverage (due to the usage of freesurfer_mapicosahedron_xhem in
    % the creation of the Caret surfaces)
    M         = caret_vol2surf_own(white.data,pial.data,images,'ignore_zeros',1);
    caret_save(fullfile(caretSDir,outfile),M); 
    fprintf('# of NaNs for %s hemi %d = %i \n',sub,h,sum(isnan(M.data)))
    M.data(isnan(M.data))=mean(M.data,'omitnan'); 
    curv_fname = [fname(1:end-4) '_' hem{h} '.curv'];
    write_curv(fullfile(outputDir,curv_fname) ,M.data,size(faces,1));
    % smooth on surface. Note that the hemisphere here
    % is lh for both hemispheres, as the indexing of
    % the rh is the same as the lh due to the
    % xhemi registration performed previously.
    % There is also a question of whether to use
    % fsaverage or fsaverage_sym as the subject here -
    % these two options are the same except for the
    % --cortex mask (medial wall) - the fsaverage_sym
    % mask is too restricted and takes out most of the
    % entorhinal cortex - so I use fsaverage!
    surf_fname_smoothed = [fname(1:end-4) '_smth' int2str(fwhm) '_' hem{h} '.mgh'];
    system(['DYLD_LIBRARY_PATH=/opt/fmrib/FreeSurfer_releases/6.0/lib/gcc/lib '...
        'mri_surf2surf --hemi lh --s fsaverage --sval ' fullfile(outputDir,curv_fname) ...
        ' --tval ' fullfile(outputDir,surf_fname_smoothed) ' --fwhm ' int2str(fwhm)]);
    
end