function O=caret_metricpermute(P,varargin)
% function caret_metricpermute(Infiles,varargin)
% Takes N files with P columns each
% and creates P files with N columns
% INPUTS: 
%   Infiles: cell array of input files, prompted if empty
% VARARGIN_OPTIONS: 
%   'inputcol',column: Only uses the column x from the input files (DEFAULT: do all)
%   'outfilenames',{'names','name'}: Name the outputfiles (DEFAULT: columnames from first metric file)
%   'outcolnames',{'names1','name2'}: Name of columns in output file(DEFAULT: names of input files) 
%   'replaceNaNs': NaNs will be replaced with 0 
% v.1.0 Joern Diedrichsen 03/06
% v.1.1 column options 
% v.2.0 : naming convention changed 
% jdiedric@bme.jhu.edu

inputcol=[];
outcolnames={}; 
outfilenames={}; 
replaceNaNs=0; 
vararginoptions(varargin,{'inputcol','outcolnames','outfilenames','replaceNaNs'},{});

if (nargin<1 | isempty(P))
    P=spm_get(inf,{'*metric'},{'Choose metric files to permute'});
end;
if (iscell(P))
    P=char(P);
end;

[a,b,type]=fileparts(P(1,:));
for i=1:size(P,1)
    filename=deblank(P(i,:)); 
    if ~exist(filename,'file') 
        warning(sprintf('File %s does not exist: replaceing with NaNs',filename)); 
        num_cols(i)=NaN;
        num_rows(i)=NaN;
    else
        M(i)=caret_load(filename);
        num_cols(i)=M(i).num_cols;
        num_rows(i)=size(M(i).data,1); 
    end;
end;

% Check if metric files are same format 
if var(num_cols(~isnan(num_cols)))>0
    warning('Number of columns is not the same');
end;
if var(num_rows(~isnan(num_rows)))>0
    error('Number of vertices are not the same');
end;

% Determine number of columns 
if (isempty(inputcol))
    inputcol=[1:min(num_cols)]; 
end; 

numrows=nanmean(num_rows); 

% Reorder into new metric files 
for file=1:length(inputcol) 
    O.num_cols=size(P,1);
    O.encoding={'BINARY'};
    O.data=zeros(numrows,length(M))*NaN; 
    for col=1:length(M)
        if (~isempty(M(col).data))           % data is present 
            O.data(:,col)=M(col).data(:,inputcol(file));
            if (~isempty(outcolnames))
                O.column_name{col}=sprintf(outcolnames{col});
            else
                [a,b,c]=fileparts(P(col,:));
                O.column_name{col}=b;
            end;
            O.column_color_mapping(col,1:2)=M(col).column_color_mapping(inputcol(file),1:2);
            if (replaceNaNs)
                O.data(isnan(O.data(:,col)),col)=0; 
            end; 
        else
            O.column_name{col}='void';
            O.column_color_mapping(col,1:2)=[0 0];
        end;
    end;    
    if (isempty(outfilenames) || length(outfilenames)<file) 
        outfilenames{file}=deblank(M(1).column_name{inputcol(file)});
        outfilenames{file}=[outfilenames{file} '.metric'];
    end;
    caret_save(outfilenames{file},O);
end;
