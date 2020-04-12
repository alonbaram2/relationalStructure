function S=caret_importfreesurfer_propatlas
% 
% Translate Freesurfer annotation file to a paint file for caret 
% Needs to translate the color structure to an areacolor file. 
% I guess you can use xmlread/xmlwrite for this 
% ---------------------------
% v1.0 Joern Diedrichsen (j.diedrichsen@ucl.ac.uk
% 

% dir='/Applications/freesurfer/subjects/fsaverage/label';
dir='/usr/local/freesurfer/subjects/fsaverage/label';

hem={'lh','rh'}; 
regions={'BA1','BA2','BA3a','BA3b','BA44','BA45','BA4a','BA4p','BA6','MT','Medial_wall','V1','V2'}; 

for h=1:2 
    D=zeros(163842,length(regions)); 
    for r=1:length(regions)
        l=read_label([],fullfile(dir,[hem{h} '.' regions{r} '.label'])); 
        D(l(:,1)+1,r)=l(:,5); 
    end; 
    M=caret_struct('metric','data',D,'column_name',regions); 
    caret_save([hem{h} '.propatlas.metric'],M); 
end;


