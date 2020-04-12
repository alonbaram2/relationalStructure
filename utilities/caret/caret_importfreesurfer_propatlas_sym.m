function S=caret_importfreesurfer_propatlas_sym
% 
% Translates the probabilisitc atlas to the symmetric 
% fsaverage_sym atlas. 
% 
% ---------------------------
% v1.0 Joern Diedrichsen (j.diedrichsen@ucl.ac.uk
% 

% dir='/Applications/freesurfer/subjects/fsaverage/label';
% dir='/usr/local/freesurfer/subjects/fsaverage/label';

sdir='/Users/jdiedrichsen/Projects/sequence_tDCS_fMRI/surfaceFreesurfer';
srcDir='/Users/jdiedrichsen/Projects/sequence_tDCS_fMRI/surfaceFreesurfer/fsaverage/surf';
targDir='/Users/jdiedrichsen/Projects/sequence_tDCS_fMRI/surfaceCaret/fsaverage_sym/LeftHem';
setenv('SUBJECTS_DIR',sdir); 

hem={'lh','rh'}; 
regions={'BA1','BA2','BA3a','BA3b','BA44','BA45','BA4a','BA4p','BA6','MT','Medial_wall','V1','V2'}; 
for r=1:length(regions)
    name{r}=sprintf('lh.%s.label',regions{r});
end; 

% Need to change source dir in routine form surf to label to make it run

% sfreesurfer_mapicosahedron_xhem('fsaverage',sdir,'curv_files',name,'surf_files',{},'hemisphere',1);

D=zeros(163842,length(regions)); 
for r=1:length(regions) 
    L=zeros(163842,1); 
    l=read_label('fsaverage',['lh.' regions{r}]); 
    L(l(:,1)+1,1)=l(:,5); 
    write_curv([srcDir filesep 'lh.' regions{r} '.curv'],L,10000);
    system(['mri_surf2surf --srcsubject fsaverage --hemi lh' ...
                    ' --trgsubject ico  --trgicoorder 7'...
                    ' --surfreg fsaverage_sym.sphere.reg'...
                    ' --srcsurfval ' 'lh.' regions{r} '.curv --src_type curv' ...
                    ' --trgsurfval ' 'lhx.' regions{r} '.curv --trg_type curv --mapmethod nnf --nsmooth-out 1']);
    D(:,r)=read_curv(fullfile(srcDir,['lhx.' regions{r} '.curv']));
end;
M=caret_struct('metric','data',D,'column_name',regions); 
caret_save(fullfile(targDir,'lh.propatlas.metric'),M); 


