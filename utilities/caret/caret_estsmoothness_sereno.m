function [FWHM]=caret_estsmoothness_sereno(cSPM,S)
% function [FWHM,resel]=caret_estsmoothness(cSPM,S)
% Estimates the smoothness on a surface 
% Using the homgenous isotropic assumtpion as used 
% reference :
% Hagler, D.J., Jr., Saygin, A.P., Sereno, M.I., 2006. 
% Smoothing and cluster thresholding for cortical surface-based group analysis of fMRI data. 
% Neuroimage 33, 1093-1103.
% ---------------------------------------------------
% v.1.0 Joern Diedrichsen j.diedrichsen@ucl.ac.uk

% ---------------------------------------------------
% remove data columns without any variance 
X=S.Nodes.data;
a=find(var(X)>0);
X=X(:,a);
[N,D]=size(X);

% ---------------------------------------------------
% calculate the residuals of the linear model and 
% standardize them 
Z=cSPM.data;
[N,n]=size(Z);
Z=(Z-cSPM.b*cSPM.X');
% STD=sqrt(nansum(Z.^2,2)/n);
% Z=Z./repmat(STD,1,n);

if (~isfield(S.Tiles,'good'))
    S.Tiles.good=true(S.num_tiles,1); 
end;

% ---------------------------------------------------
% Generate all possible M pairs of nodes 
Z1=Z(S.Edges.data(:,1),:); 
Z2=Z(S.Edges.data(:,2),:); 

Zd=[Z1-Z2];
Xd=[S.Edges.Length]; 

j=find(~isnan(sum(Zd,2))); 
Zd=[Zd(j,:);-Zd(j,:)]; 
Xd=[Xd(j,:);Xd(j,:)]; 
Z=[Z1(j,:);Z2(j,:)];

FWHM=mean(Xd)*sqrt(-2*log(2)/log(1-var(Zd(:))/(2*var(Z(:)))));
