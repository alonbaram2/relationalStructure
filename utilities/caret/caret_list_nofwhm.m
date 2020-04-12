function [TabDat,threshold_con]=caret_list_nofwhm(S,cSPM,u,k,varargin);
% function TabDat=caret_list_nofwhm(Surface,cSPM,u,k,varargin);
% Function that generates a Table with clusters, local maxima
% and corrected P-values for a certain threshold u, k 
% Usage:
% Tabdat=caret_list(Surface,cSPM,u,k,varargin)
%       Surface: a surface structure (see caret_getsurface)
%       cSPM: statistical map structure (see caret_getcSPM)
%       u: height-threshold
%       k: size-threshold (mm)
%       varargin:
%           'sort_by',{'p_max'/'area'/'p_cluster'/'x_coord'}
%           'contrast',contrast-number (default is 1)
%           'sign',                     : Signed t-test +1 / -1
%           'coord',avrg_PALS_file      : Looks up the coordinates on a different
%                                       Surface
%           'mask',vec                  : logical mask that defines a subregion of the surface
%                                         that is used as a search region
%           'Nmask',prop                : mask for proportion of subjects for which you need observations (default: 0.7)
%           'save_thresholded',fname,col: : Saves the resulting masked map
%                                         as a column to a new or exisiting
%                                         metric file
% caret_list('txtlist',TabDat)
%       generates a text print-out of the table
%
% OUTPUT: The table has the following columns:
% c	NumN Area max(Z) p(unc) p(cor) p(cl.) X(mm) Y(mm) Z(mm) Node
%   c: number of cluster
%   NumN: Number of Nodes in cluster
%   Area: area of cluster in mm
%   max(Z): maximum of statistical value in that cluster
%   p(unc): uncorrected p-value for maximum
%   p(cor): corrected p-value for maximum
%   p(cl.): cluster-wise p-value for the cluster of that size, based on
%           permutation stats 
%   X,Y,Z:  X-,Y-, and Z- coordinate of the maximum in the space the map is
%           in. Only if the map is based on a fiducial in a MNI-space do these
%           correspond to MNI-coordinates
%   Node :  Number of Node at Maximum 
%           This is 1-based, so you need to subtract 1 if you use caret to identify 
%           the node over Surface->Identify
%-----------------------------------------------------------------------
% v.2.0. Joern Diedrichsen 12/4/10 j.diedrichsen@ucl.ac.uk

sort_by='p_max';
contrast=1;
sign=1;
coord=[];
mask=[];        % binary masking image for small subsets
Nmask=0.7;      % Automatically masking for number of data avaiable
save_thresholded=[]; 

vararginoptions(varargin,{'sort_by','contrast','sign','coord','mask','save_thresholded','Nmask'});

% Mask the image and adjust the search region
con=cSPM.con(contrast);
if (sign==-1)
    con.Z=-con.Z;
end;


if (isempty(mask))
    mask=cSPM.N>size(cSPM.data,2)*Nmask;
else
    mask=(cSPM.N>size(cSPM.data,2)*Nmask) & mask;
end;

if (~isempty(coord))
    C=caret_load(coord);
end;

%-----------------------------------------------------------------------
% Calculate Clusters and local maxima
A=caret_clusters(S,con.Z>u & mask);
a=unique(A);
cluster=0;
TabDat.dat=[];
threshold_con=zeros(size(con.Z));
for j=1:length(a)-1
    i=find(A==j);
    area=sum(S.Nodes.area(i));
    if (area./(fwhm^2))>k
        threshold_con(i)=sign*con.Z(i);  % If cluster is big enough- store thresholded contrast
        cluster=cluster+1;
        TabDat.dat(cluster,1)=cluster;                  % Number of cluster
        TabDat.dat(cluster,2)=length(i);                % Number of nodes in cluster
        TabDat.dat(cluster,3)=area;                     % Total area of cluster
        [maxZ,maxInd]=max(con.Z(i));   % Maximal T-value of cluster
        TabDat.dat(cluster,4)=sign*maxZ;
        TabDat.dat(cluster,5)=caret_P(1,0,maxZ,con.df,con.STAT,1); % Uncorrected p-vale
        TabDat.dat(cluster,6)=caret_P(1,0,maxZ,con.df,con.STAT,R); % Corrected p-value
        TabDat.dat(cluster,7)=caret_P(1,TabDat.dat(cluster,3)/(fwhm^2),u,con.df,con.STAT,R);
        if (~isempty(coord))
            TabDat.dat(cluster,8:10)=C.data(i(maxInd),:);
            TabDat.dat(cluster,11)=i(maxInd);
        else
            TabDat.dat(cluster,8:10)=S.Nodes.data(i(maxInd),:);
            TabDat.dat(cluster,11)=i(maxInd);
        end;
    end;
end;
Pz              = caret_P(1,0,u,con.df,con.STAT,1);
Pu              = caret_P(1,0,u,con.df,con.STAT,R);
[P Pn Em En EN] = caret_P(1,k,u,con.df,con.STAT,R);

%[P,p,Em,En,EN] = caret_P(num_clusters,4,u,cSPM.df,cSPM.STAT,R);
%-----------------------------------------------------------------------
%-Headers for text table...
TabDat.tit = cSPM.title;
TabDat.hdr = {'  c','NumN','Area','max(Z)','p(unc)','p(cor)','p(cl.)','X(mm)','Y(mm)','Z(mm)','Node'};
TabDat.fmt = {'%3d','%4d','%-3.2f','%0.3f','%0.3f', '%0.3f', '%0.3f','%4.2f','%4.2f','%4.2f','%6d'};				%-XYZ

%-Footnote with SPM parameters
%-----------------------------------------------------------------------
TabDat.ftr    = cell(4,2);
TabDat.ftr{1} = ...
    sprintf('Height threshold: %c = %0.2f, p = %0.3f (%0.3f)',...
    con.STAT,u,Pz,Pu);
TabDat.ftr{2} = ...
    sprintf('Extent threshold: k = %0.2f sqmm, p = %0.3f (%0.3f)',...
    k*fwhm^2,Pn,P);
TabDat.ftr{3} = ...
    sprintf('Expected sqmm per cluster, <k> = %0.3f',En*fwhm.^2);
TabDat.ftr{4} = ...
    sprintf('Expected number of clusters, <c> = %0.2f',Em*Pn);
TabDat.ftr{5} = ...
    sprintf('Expected surface above threshold (sqmm) = %0.2f',Em*En*fwhm.^2);
TabDat.ftr{6} = ...
    sprintf('Degrees of freedom = [%0.1f, %0.1f]',con.df);
TabDat.ftr{7} = ...
    sprintf('Smoothness FWHM = %0.1f (mm) ',fwhm);
TabDat.ftr{8} = ...
    sprintf('Search vol: EC:%d Circ:%3.0fmm Surf:%4.0fmm2; %0.1f resels',Area(1),Area(2),Area(3),R(end));

% -----------------------------------------------------------------
% Put all the necessary information in Table-strcuture
TabDat.Em=Em;
TabDat.En=En*fwhm.^2;
TabDat.Pn=Pn;
TabDat.EN=EN*fwhm.^2;
if (isempty(TabDat.dat))
    TabDat.TotA=0;
    TabDat.m=0;
    TabDat.n=0;
else
    TabDat.TotA=sum(TabDat.dat(:,3));
    TabDat.m=cluster;
    TabDat.n=mean(TabDat.dat(:,3));
    % -----------------------------------------------------------------
    % Sort and list the Table
    switch (sort_by)
        case 'p_max'
            sort_col=-TabDat.dat(:,4);
        case 'area'
            sort_col=-TabDat.dat(:,3);
        case 'p_cluster'
            sort_col=TabDat.dat(:,7);
        case 'x_coord'
            sort_col=TabDat.dat(:,8);
    end;
    [A,IND]=sort(sort_col);
    TabDat.dat=TabDat.dat(IND,:);
    caret_listprint('txtlist',TabDat);
end;

if (~isempty(save_thresholded)) 
    [direct,fname,ext,col]=spm_fileparts(save_thresholded); 
    col=str2num(col(2:end)); 
    fname=fullfile(direct,[fname ext]); 
    if (exist(fname,'file')) 
        D=caret_load(fname); 
        D.column_name{col}=['thresholded']; 
        D.num_cols=max(D.num_cols,col); 
        D.data(:,col)=threshold_con;
        D.column_color_mapping(col,1:2)=[u max(threshold_con)]; 
        caret_save(fname,D); 
    else 
        error('not implemented yet'); 
    end; 
end; 
