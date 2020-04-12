function [TabDat,threshold_con]=caret_list(S,cSPM,u,k,varargin);
% function TabDat=caret_list(Surface,cSPM,u,k,varargin);
% Function that generates a Table with clusters, local maxima
% and corrected P-values for a certain threshold u
% Usage:
% Tabdat=caret_list(Surface,cSPM,u,k,varargin)
%       Surface: a surface structure (see caret_getsurface)
%       cSPM: statistical map structure (see caret_getcSPM)
%       u: height-threshold (if u<1 it is interpreted as uncorrected p) 
%       k: size-threshold(mm2) if k<1 it is interpreted as corrected p cluster) 
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
%           'label',paintname           : Name of paint file, whos entries
%                                         are used to label the 
%           'Area', [EC length area]    : Area to correct over: use this 
%                                         to correct for both hemispheres 
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
%   p(cl.): cluster-wise p-value for the cluster of that size
%   X,Y,Z:  X-,Y-, and Z- coordinate of the maximum in the space the map is
%           in. Only if the map is based on a fiducial in a MNI-space do these
%           correspond to MNI-coordinates
%   Node :  Number of Node at Maximum 
%           This is 1-based, so you need to subtract 1 if you use caret to identify 
%           the node over Surface->Identify
%           if label paint file is given, the name of the paint is listed 
%-----------------------------------------------------------------------
% v.1.0 Joern Diedrichsen 12/10/04 jdiedric@bme.jhu.edu
% v.2.0. 12/4/10 j.diedrichsen@ucl.ac.uk
% v.3.0 2/2/11: Uses now fMRIstats' stat_threshold routine for better
%               P-values 

c=5;
sort_by='p_max';
contrast=1;
sign=1;
coord=[];
mask=[];        % binary masking image for small subsets
Nmask=0.7;      % Automatically masking for number of data 
Pcorr=0.05; 
save_thresholded=[]; 
label=[]; 
Area=[]; 

vararginoptions(varargin,{'sort_by','contrast','sign','coord','mask',...
    'save_thresholded','Nmask','label','Area'});

% Load label where applicable
if (ischar(label))
    label=caret_load(label); 
end; 

% If not done yet, calculate the average areas of the Tiles and Nodes
if ~isfield(S,'A');
    fprintf('Calculate Search Area\n');
    S=caret_calcarea(S);
end;

% Mask the image and adjust the search region
con=cSPM.con(contrast);
if (isempty(mask))
    mask=cSPM.N>size(cSPM.data,2)*Nmask;
else
    mask=(cSPM.N>size(cSPM.data,2)*Nmask) & mask;
end;
if (isempty(Area))
    Area=caret_calc_EC(S,mask);
end; 

% Estimate the smoothness of the map
if ~isfield(cSPM,'fwhm');
    fprintf('Calculate Smoothness\n');
    fwhm=caret_estsmoothness_sereno(cSPM,S); 
else 
    fwhm=cSPM.fwhm; 
end; 

% If negative contrast, turn statisitics around (only for T) 
if (sign==-1)
    con.Z=-con.Z;
end;

% If different coordinate file is given, load that 
if (~isempty(coord))
    C=caret_load(coord);
end; 

% Get corrected thresholds 
switch (con.STAT)
    case 'T' 
        df=[con.df(2) 0];
        if (u<1) 
            Punc=u;
            u=tinv(1-Punc,con.df(2)); 
        else 
            Punc=1-tcdf(u,con.df(2));
        end;         
    otherwise 
        error('not implemented yet'); 
end; 

% At the uncorrected threshold, what would be the critical values
if (k<1)   
    Pcl=k;
    [Ucrit,k]=stat_threshold(Area,sum(mask),fwhm,df,Pcorr,u,Pcl);
else 
    [Ucrit,Pcl]=stat_threshold(Area,sum(mask),fwhm,df,Pcorr,u,k);    
end;
[Ucrit,Kcrit]=stat_threshold(Area,sum(mask),fwhm,df,0.05,u,0.05);
[P Pn Em En EN] = caret_P(1,k/(fwhm.^2),u,con.df,con.STAT,Area./(fwhm.^[0 1 2]));

%-----------------------------------------------------------------------
% Calculate Clusters and local maxima
A=caret_clusters(S,con.Z>u & mask);
a=unique(A);
cluster=0;
TabDat.dat={};
threshold_con=zeros(size(con.Z));
for j=1:length(a)-1
    i=find(A==j);
    area(j)=sum(S.Nodes.area(i));
    if area(j)>k
        threshold_con(i)=sign*con.Z(i);  % If cluster is big enough- store thresholded contrast
        cluster=cluster+1;
        TabDat.dat{cluster,1}=cluster;                  % Number of cluster
        TabDat.dat{cluster,2}=area(j);                     % Total area of cluster
        [maxZ,maxInd]=max(con.Z(i));                    % Maximal T-value of cluster
        TabDat.dat{cluster,3}=sign*maxZ;
        
        % Now call fMRI stat to get corrected thresholds 
        [pcorr,pcl]=stat_threshold(Area,sum(mask),fwhm,df,maxZ,u,area(j));
        punc=1-tcdf(maxZ,con.df(2));
        TabDat.dat{cluster,4}=punc; % Uncorrected p-vale
        TabDat.dat{cluster,5}=pcorr; % Corrected p-value
        TabDat.dat{cluster,6}=pcl;  % Cluster-wise threshold 
        if (~isempty(coord))
            TabDat.dat{cluster,7}=C.data(i(maxInd),1);
            TabDat.dat{cluster,8}=C.data(i(maxInd),2);
            TabDat.dat{cluster,9}=C.data(i(maxInd),3);
        else
            TabDat.dat{cluster,7}=S.Nodes.data(i(maxInd),1);
            TabDat.dat{cluster,8}=S.Nodes.data(i(maxInd),2);
            TabDat.dat{cluster,9}=S.Nodes.data(i(maxInd),3);
        end;
        if (isempty(label))
            TabDat.dat{cluster,10}=i(maxInd);
        else 
            TabDat.dat{cluster,10}=label.paintnames{label.data(i(maxInd))};
        end; 
    end;
end;



%[P,p,Em,En,EN] = caret_P(num_clusters,4,u,cSPM.df,cSPM.STAT,R);
%-----------------------------------------------------------------------
%-Headers for text table...
TabDat.tit = cSPM.title;
if (isempty(label))
    TabDat.hdr = {'  c','Area','max(Z)','p(unc)','p(cor)','p(cl.)','X(mm)','Y(mm)','Z(mm)','Node'};
    TabDat.fmt = {'%3d','%-3.2f','%-0.2f','%0.3f', '%0.3f', '%0.3f','%4.2f','%4.2f','%4.2f','%6d'};				%-XYZ
else 
    TabDat.hdr = {'  c','Area','max(Z)','p(unc)','p(cor)','p(cl.)','X(mm)','Y(mm)','Z(mm)','Region'};
    TabDat.fmt = {'%3d','%-3.2f','%-0.2f','%0.3f', '%0.3f', '%0.3f','%4.2f','%4.2f','%4.2f','%s'};				%-XYZ
end;    

%-Footnote with SPM parameters
%-----------------------------------------------------------------------
TabDat.ftr    = cell(4,2);
TabDat.ftr{1} = ...
    sprintf('Height threshold: %c = %0.2f, p = %0.3f ',...
    con.STAT,sign*u,Punc);
TabDat.ftr{2} = ...
    sprintf('Height threshold crit: %c = %0.2f, P = %0.3f',...
    con.STAT,sign*Ucrit,0.05);
TabDat.ftr{3} = ...
    sprintf('Extent threshold: k = %0.2f sqmm, P = %0.3f',...
    k,Pcl);
TabDat.ftr{4} = ...
    sprintf('Extent threshold crit: k = %0.2f sqmm P = 0.05',...
    Kcrit);
TabDat.ftr{5} = ...
    sprintf('Expected sqmm per cluster, <k> = %0.3f',En*fwhm.^2);
TabDat.ftr{6} = ...
    sprintf('Expected number of clusters, <c> = %0.2f',Em);     % There was a pn here which was likely not correct
TabDat.ftr{7} = ...
    sprintf('Expected surface above threshold (sqmm) = %0.2f',Em*En*fwhm.^2);
TabDat.ftr{8} = ...
    sprintf('Degrees of freedom = [%0.1f, %0.1f]',con.df);
TabDat.ftr{9} = ...
    sprintf('Smoothness FWHM = %0.1f (mm) ',fwhm);
TabDat.ftr{10} = ...
    sprintf('Search vol: EC:%d Circ:%3.0fmm Surf:%4.0fmm2; %0.1f resels',Area(1),Area(2),Area(3),Area(3)./(fwhm.^2));

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
    caret_listprint('txtlist',TabDat);
else
    TabDat.TotA=sum(area);
    TabDat.m=cluster;
    TabDat.n=mean(area);
    % -----------------------------------------------------------------
    % Sort and list the Table
    switch (sort_by)
        case 'p_max'
            sort_col=-vertcat(TabDat.dat{:,3});
        case 'area'
            sort_col=-vertcat(TabDat.dat{:,2});
        case 'p_cluster'
            sort_col=vertcat(TabDat.dat{:,7});
        case 'x_coord'
            sort_col=vertcat(TabDat.dat{:,8});
    end;
    [A,IND]=sort(sort_col);
    for i=1:length(IND)
        for j=1:size(TabDat.dat,2) 
            newdat{i,j}=TabDat.dat{IND(i),j};
        end; 
    end; 
    TabDat.dat=newdat;
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
