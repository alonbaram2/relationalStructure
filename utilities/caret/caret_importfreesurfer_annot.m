function [P,AC]=caret_importfreesurfer_annot(annot_name,paint_name,areacolor);
% 
% Translate Freesurfer annotation file to a paint file for caret 
% Needs to translate the color structure to an areacolor file. 
% I guess you can use xmlread/xmlwrite for this - looks very involved, though!  
% ---------------------------
% v1.0 Joern Diedrichsen (j.diedrichsen@ucl.ac.uk
% 
[v,l,ct]=read_annotation(annot_name); 
COLOR=ct.table(:,1:3);            % These are the color values 
VALUES=ct.table(:,5); 
% Relabel them from 1-36 
for i=1:length(VALUES)
    l(l==VALUES(i))=i;
end; 

P=caret_struct('paint','data',l,'paintnames',ct.struct_names);
caret_save(paint_name,P); 


% Make the area-color file (comma-seperated) 
AC{1}='CSVF-FILE,0,,,,,,,';
AC{2}='csvf-section-start,header,2,,,,,,';
AC{3}='tag,value,,,,,,,';
AC{4}='Caret-Version,5.65,,,,,,,';
AC{5}='Date,2012-08-19T11:19:34,,,,,,,';
AC{6}='comment,,,,,,,,';
AC{7}='encoding,COMMA_SEPARATED_VALUE_FILE,,,,,,,';
AC{8}='pubmed_id,,,,,,,,';
AC{9}='csvf-section-end,header,,,,,,,';
AC{10}='csvf-section-start,Colors,9,,,,,,';
AC{11}='Name,Red,Green,Blue,Alpha,Point-Size,Line-Size,Symbol,SuMSColorID';
for i=1:size(COLOR,1)
    AC{end+1}=sprintf('%s,%d,%d,%d,255,2.000,1.000,POINT,',ct.struct_names{i},COLOR(i,1),COLOR(i,2),COLOR(i,3));
end; 
AC{end+1}='csvf-section-end,Colors,,,,,,,'; 

fid=fopen(areacolor,'w+'); 
for i=1:length(AC); 
    fprintf(fid,'%s\n',AC{i}); 
end; 
fclose(fid); 