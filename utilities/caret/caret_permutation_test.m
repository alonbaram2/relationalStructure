function [u,k,U,K]=caret_permutation_test(S,cSPM,u,varargin);
% function [u,k,U,K]=caret_permutation_test(S,cSPM,u,varargin);
% caret_list('txtlist',TabDat)
%       generates a text print-out of the table
%
% OUTPUT: The table has the following columns:
%   u: height threshold at p<0.05
%   k: Cluster size area (mm2) at p<0.05
%   U: List of maxima 
%   K: List of biggest cluster 
%-----------------------------------------------------------------------
% v.2.0. Joern Diedrichsen 12/4/10 j.diedrichsen@ucl.ac.uk

contrast=1;
mask=[];        % binary masking image for small subsets
Nmask=0.9;      % Automatically masking for number of data avaiable

vararginoptions(varargin,{'contrast','sign','mask','Nmask'});

% Mask the image and adjust the search region
con=cSPM.con(contrast);

if (isempty(mask))
    mask=cSPM.N>size(cSPM.data,2)*Nmask;
else
    mask=(cSPM.N>size(cSPM.data,2)*Nmask) & mask;
end;

%-----------------------------------------------------------------------
% Find all permutations 
N=size(cSPM.data,2); 
if (N>13)
    error('completet permutation test takes too long'); 
end; 
X=ones(2^(N-1),N);
M=2^(N-1);
for n=1:N-1
    A=[ones(M/(2^n),1);-ones(M/(2^n),1)];
    X(:,n+1)=repmat(A,2^(n-1),1);
end; 

for m=1:M 
    D=bsxfun(@times,cSPM.data,X(m,:)); 
    b=nanmean(D,2);
    ResVar=nanvar(D,0,2);
    ResVar(ResVar==0)=NaN;
    Z=b./sqrt(cSPM.ResVar/N);  
    A=caret_clusters(S,Z>u & mask);
    a=unique(A);
    for j=1:length(a)-1
        i=find(A==j);
        area(j)=sum(S.Nodes.area(i));
    end; 
    U(m,1)=max(Z(mask,1)); 
    K(m,1)=max(area); 
end; 
u=prctile(U,95); 
k=prctile(K,95); 


