function Out=caret_smooth(P,varargin);
% Smooth metric files 
% Out=caret_smooth(P,varargin);
% Uses caret_command to smooth the files 
% INPUT: 
%   P: cell array of metric files to smooth 
% VARARGIN: 
%   'coord',file:   Coordinate file 
%   'topo',topo:    Topology file 
%   'algorithm',str: Smoothing algorithm: default AN 
%           AN: Average Neighbour 
%           DILATE: Dilatiom
%           FWHM: Full-Width-Half Maximum
%           GAUSS: Gaussian volume smoothing 
%           GEOGAUSS: Geodesic Gaussian algorithm 
%           http://brainvis.wustl.edu/wiki/index.php/Caret:Documentation:MetricSmoothing
%   'gauss',sigma
%   'fwhm',FWHM:    desired FWHM for fwhm iteration
% OUTPUT: 
%   Cell array of smooth metric files 
% v1.0 Joern Diedrichsen (j.diedrichsen@ucl.ac.uk)
%
% Tobias Wiestler (31.05.2011): Force to write binary metric files 
prefix='s';
iterations=10;
strength=0.5;
algorithm='AN';
fwhm=6; 
coord=[];
topo=[];
vararginoptions(varargin,{'prefix','iterations','strength','algorithm','coord','topo','fwhm'});
if (nargin<1 | isempty(P))
    P=spm_get(inf,'*.metric','Pick Metric files to smooth');
elseif (iscell(P))
    P=char(P);
end;
if (isempty(coord))
    coord=spm_get(1,'*.coord','Pick Coordinate file');
end;
if (isempty(topo))
    topo=spm_get(1,'*.topo','Pick Topo file');
end;
N=size(P,1);
for i=1:N 
    [dir,name,postfix]=fileparts(deblank(P(i,:)));
    Out{i}=fullfile(dir,[prefix name postfix]); 
    switch (algorithm) 
        case 'AN'
            comm=sprintf('caret_command -metric-smoothing %s %s %s %s %s %d %d -WRITE-FILE-FORMAT-METRIC BINARY',...
            coord,topo,fullfile(dir,[name postfix]),Out{i},algorithm,iterations,strength);
        case 'FWHM'
            comm=sprintf('caret_command -metric-smoothing %s %s %s %s %s %d %d -fwhm %f -WRITE-FILE-FORMAT-METRIC BINARY',...
            coord,topo,fullfile(dir,[name postfix]),Out{i},algorithm,1000,1,fwhm);
        otherwise 
            error ('not implememented yet'); 
    end; 
    fprintf('%s\n',comm) 
    [err,out]=system(comm);
    fprintf('Smoothed file %s\n%d:%s\n',[name postfix],err,out);
end;

