function [PN,AC]=caret_combinePaint(paint_name,outpaintname,outareacolor,varargin);
% Combines Paint ROIs to a new ROIS 
% The default makes Lobular ROIs from the desikan freesurfer atlas 

P=caret_load(paint_name); 

% Name of new ROIs 
names={'frontal','parietal','occipital','temporal'};

% Areas included in each ROI 
areas{1}={'precentral','caudalmiddlefrontal','parsopercularis','superiorfrontal',...
    'parstriangularis','paracentral'...
    'posteriorcingulate','caudalanteriorcingulate','rostralmiddlefrontal',...
    'caudalanteriorcingulate','rostralanteriorcingulate',...
    'medialorbitofrontal','lateralorbitofrontal','parsorbitalis'}; 
areas{2}={'postcentral','superiorparietal','supramarginal','inferiorparietal','superiorparietal','precuneus','isthmuscingulate'}; 
areas{3}={'lateraloccipital','cuneus','pericalcarine','lingual'};
areas{4}={'superiortemporal','middletemporal','inferiortemporal',...
    'bankssts','transversetemporal','insula','fusiform',...
    'parahippocampal','entorhinal','temporalpole'}; 

% COLORs for new ROIS 
colors=[255 0 0;0 255 0;0 0 255;160 160 160]; 

vararginoptions(varargin,{'areas','colors','names'}); 

% initialize the data field 
Data=P.data*0; 

% Assign indices 
for i=1:length(areas) 
    for j=1:length(areas{i})
        indx=find(strcmp(areas{i}{j},P.paintnames)); 
        if (isempty(indx))
            error(sprintf('%s not found in paint file',areas{i}{j})); 
        end; 
        Data(P.data==indx)=i; 
    end; 
end; 

PN=caret_struct('paint','data',Data,'paintnames',names'); 
caret_save(outpaintname,PN);

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
for i=1:length(areas)
    AC{end+1}=sprintf('%s,%d,%d,%d,255,2.000,1.000,POINT,',names{i},colors(i,1),colors(i,2),colors(i,3));
end; 
AC{end+1}='csvf-section-end,Colors,,,,,,,'; 

fid=fopen(outareacolor,'w+'); 
for i=1:length(AC); 
    fprintf(fid,'%s\n',AC{i}); 
end; 
fclose(fid); 