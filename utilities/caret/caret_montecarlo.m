function [E,cSPM,S]=caret_montecarlo(S,u,N,fwhm,varargin)
% function X=caret_montecarlo(S,cSPM)
% caret_montecarlo(S,u,fwhm);   generates the necessary smoothed randomization
%                            files on the surface
% Monte-carlo study on surface statistics

num_iter=1; 
mask=[]; 
check_caretlist=1; 
vararginoptions(varargin,{'num_iter','mask'});

if (isempty(mask))
    mask=true(S.num_nodes,1); 
end; 

[dir]=fileparts(which('caret_montecarlo'));
surf_dir=fullfile(dir,'Surfaces');
fname=sprintf('random_data_%2.1f.metric',fwhm);

S=caret_calcarea(S);


if (~exist(fullfile(surf_dir,['s' fname])))
    M=caret_struct('metric','data',normrnd(0,1,S.num_nodes,60));
    caret_save(fname,M);
    caret_smooth({fname},'coord',S.coord_file,'topo',S.topo_file,'algorithm','FWHM','fwhm',fwhm);
end;

M=caret_load(fullfile(surf_dir,['s' fname])); 
P=size(M.data,2); 
for i=1:num_iter
    indx=sample_wr([1:P],N);
    D=M.data(:,indx); 
    x=unidrnd(2,1,N)*2-3; 
    D=bsxfun(@times,D,x); 
    
    % Do a one-sample t-test on this 
    %b=nanmean(D,2);
    %ResVar=nanvar(D,0,2);
    %ResVar(ResVar==0)=NaN;
    %Z=b./sqrt(ResVar/N);  
    cSPM=caret_getcSPM('onesample_t','data',D); 
    Z=cSPM.con.Z;
    
    A=caret_clusters(S,Z>u & mask);
    
    % Find clusters 
    a=unique(A);
    area=zeros(length(a),1);
    for j=1:length(a)-1
        in=find(A==j);
        area(j)=sum(S.Nodes.area(in));
    end; 
    E.U(i,1)=max(Z(mask,1)); 
    E.K(i,1)=max(area); 
    E.u(i,1)=u; 
    E.N(i,1)=N; 
    if (check_caretlist)
        E.numcluster(i,1)=length(a)-1;
        E.totalarea(i,1)=sum(area);
        E.numnodes(i,1)=sum(Z>u & mask); 
        FWHM=caret_estsmoothness(cSPM,S)'; 
        E.fwhm(i,1)=(prod(FWHM).^(1./3));
        E.fwhm_sereno(i,1)=caret_estsmoothness_sereno(cSPM,S);
    end; 
    
    % Tab=caret_list(S,cSPM,u,10,'sort_by','area');
    fprintf('.'); 
end;
fprintf('\n');
umax=prctile(E.U,95); 
k=prctile(E.K,95); 
