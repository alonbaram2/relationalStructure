function Y=caret_crosssection(border,datafile,varargin);
% Samples a cross-section from a metric file 
% INPUT: 
%   border: Borderproj file:
%   Datafiles: name of metric, coord, paint or RGB file to be sampled 
%           Data to be maps
% VARARGIN: 
%   size: radius around each node (default = 0) 
% OUTPUT: 
%   Y:   data 
% To get the coordinates for the points, simply sample from the respective
% coord file. 
% (C) Joern Diedrichsen 2013
[a,b,c,d]=spm_fileparts(border);
if (~strcmp(c,'.borderproj'))
        error('need a border projection file');
end;
B=caret_load(border); 
vertex=B.Border.vertex; 

D=caret_load(datafile);

for i=1:3 
    Y(:,:,i)=D.data(vertex(:,i),:); 
end; 
Y=mean(Y,3);