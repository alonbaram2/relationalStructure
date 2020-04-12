function X=caret_testsurfacestats(varargin)
% Monte-carlo study on surface statistics

X=[];
surfacefile='lh.avrgsurface.mat';
coordfile=fullfile(['lh.FIDUCIAL.coord']);
topofile=fullfile(['lh.CLOSED.topo']);

% Smoothing
iterations=10;
strength=1;
N=10;       % number of subjects
n=20;        % Number of iterations


load(surfacefile);


for i=1:n
    M=caret_struct('metric','data',normrnd(0,1,S.num_nodes,N));
    caret_save('test.metric',M);
    
    sm=caret_smooth('test.metric','coord',coordfile,'topo',topofile,'iterations',iterations,'strength',strength);
    cSPM=caret_getcSPM('onesample_t','data',sm{1},[]);
    T=caret_list(S,cSPM,3,0);
    D.numClusterSig=sum(T.dat(:,7)<0.05);
    D.numPeaksSig=sum(T.dat(:,6)<0.05);
    D.meanArea=mean(T.dat(:,3));
    D.numClusters=length(T.dat);
    
    X=addstruct(X,D);
end;


