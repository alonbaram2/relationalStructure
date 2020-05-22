
rootData    = '/home/fs0/abaram/scratch/NN2020';
subjects = {'sub-01','sub-02','sub-03','sub-04','sub-05','sub-06','sub-07','sub-08','sub-09',...
    'sub-10','sub-11','sub-12','sub-13','sub-14','sub-15','sub-16','sub-17','sub-18','sub-19',...
    'sub-20','sub-21','sub-22','sub-23','sub-24','sub-25','sub-26','sub-27','sub-28'};

glm = 'GLM3';
fwhm = 5;

for iSub=1:length(subjects)
    sub = subjects{iSub};
    
    spmDir = fullfile(rootData,'spmDirs',glm,sub);
    outputDir = fullfile(spmDir,'surfCons');
    if ~exist(outputDir,'dir')
        mkdir(outputDir)
    end
    cd(spmDir);    
%     filenames = dir('con*.nii');
    filenames = dir('con*.nii');
    
    hemName = {'LeftHem','RightHem'};
    hem     = {'lh','rh'};
    caretDir = fullfile(rootData,'caretDirs',sub);
    [~, faces] = read_surf(fullfile(rootData,'freesurfer-subjects','fsaverage','surf','lh.sphere.reg'));
    for iBeta = 1:length(filenames)
        fname = filenames(iBeta).name;
        for h=1:2
            % the surfaces in caretSDir are already in world (not
            % RAS) coordinates
            caretSDir = fullfile(caretDir,sprintf(['x',sub]),hemName{h});
            white     = caret_load(fullfile(caretSDir,[hem{h} '.WHITE.coord']));
            pial      = caret_load(fullfile(caretSDir,[hem{h} '.PIAL.coord']));
            images    = fullfile(spmDir,fname);
            outfile   = [fname(1:end-4) '.metric'];
            % map data onto surface, according to world coords.
            % notice vertex indeces correspond (for both hemi) to the left hemi of
            % fsaverage (due to the usage of freesurfer_mapicosahedron_xhem in
            % the creation of the Caret surfaces)
            M         = caret_vol2surf_own(white.data,pial.data,images,'ignore_zeros',1);
            caret_save(fullfile(caretSDir,outfile),M);
            fprintf('# of NaNs for %s hemi %d = %i \n',sub,h,sum(isnan(M.data)))
            M.data(isnan(M.data))=mean(M.data,'omitnan');
            surf_fname = [fname(1:end-4) '_' hem{h}];
            write_curv(fullfile(outputDir,surf_fname) ,M.data,size(faces,1));
            % smooth on surface. Note that the hemisphere here
            % is lh for both hemispheres, as the indexing of
            % the rh is the same as the lh of fsaverage_sym.
            % There is also a question of whether to use
            % fsaverage or fsaverage_sym as the subject here -
            % these two options are the same except for the
            % --cortex mask (medial wall) - the fsaverage_sym
            % mask is too restricted and takes out most of the
            % entorhinal cortex - so I use fsaverage!
            surf_fname_smoothed = [fname(1:end-4) '_smth' int2str(fwhm) '_' hem{h} '.mgh'];
            system(['DYLD_LIBRARY_PATH=/opt/fmrib/FreeSurfer_releases/6.0/lib/gcc/lib '...
                'mri_surf2surf --hemi lh --s fsaverage --sval ' fullfile(outputDir,surf_fname) ...
                ' --tval ' fullfile(outputDir,surf_fname_smoothed) ' --fwhm ' int2str(fwhm)]);            
        end
    end
end