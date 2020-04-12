function COG=caret_cog(Y,c,coord,ROI,varargin); 
% function varargout=caret_cog(Y,C,coord,varargin)
% INPUT: 
%    Y       : Data, a VxQ matrix for V vertices and Q conditions*subjects
%    c       : 1xQ condition labels 
%    coord   : File name for the coord file (usually flat)
%    ROI     : Indicator matrix giving the ROI 
% OUTPUT: 
%    COG       : Qx2/3 vectore of COG coordinates
% VARARGIN: 
%    'measure'            : 'cog' (default) or 'max' 
%                           cog calculates the center of gravity (igonoring
%                           negative numbers) 
%                           max takes the location of the maximum for each
%                           subject - makes sense only on smoothed data 
%   for plotting: 
%    'topo',filename      : Topology file: if given it will plot  
%    'xlims', [xmin xmax] : xlimiters for plotting (default: the ROI)
%    'ylims', [ymin ymax] : ylimiters for plotting (default: the ROI)
%    'underlay',data      : data for color underlay
%    'uscale',[min max]   : Min (black) and max (white) of the underlay
%    'dscale'             : [threshold(min) and peak(max)] for display
%    'doffset',num        : Brightness of the smallest values 
%    'border'             : Possible adding a border
%    'alpha'              : Alpha value - blending for 
%   for testing 
%    'test',{'paired','independent'}: If given, the function conducts a
%                       Hotellings T^2 test too see if the COGs of the to groups (independnet)
%                       or conditions are different from each other 
% (C) j.diedrichsen@ucl.ac.uk, 2014 

topo=[]; 
xlims=[]; 
ylims=[]; 
underlay=[]; 
uscale=[]; 
dscale=[];
doffset=0.4; 
border=[]; 
alpha=1; 
markersize=6; 
threshold = []; 
color=[1 3 2]; % 1 red: 2: green 3:blue  
colorC={'r','g','b'}; 
plotmean = 0; 
test=''; 

vararginoptions(varargin,{'measure','topo','xlims','ylims','data','underlay','uscale',...
                          'M','dscale','border','alpha','color','markersize','test','threshold','plotmean'});
colorC={colorC{color}}; % Adjust the colors for symbol plotting 
C=caret_load(coord);
    
% Get the XYZ coordinates for stuff in the ROI 
XYZ=C.data(ROI,:); 
Q=size(Y,2); 

% Determine the COG for each of the of the columns 

% Enforce positivity for COG calculation 
switch (measure) 
    case 'cog' 
        Y(Y<0)=0; 
        for q=1:Q 
            A=bsxfun(@times,Y(ROI,q),XYZ); 
            COG(q,:)=sum(A)./sum(Y(ROI,q)); 
        end; 
    case 'max'
        for q=1:Q 
            [~,idx]=max(Y(ROI,q)); 
            COG(q,:)=XYZ(idx,:); 
        end; 
end;         

% Make a nice RGB plot if so desired 
if (~isempty(topo)) 
    MIN=min(XYZ); 
    MAX=max(XYZ); 
    if (isempty(xlims))
        xlims=[MIN(1) MAX(1)]; 
    end; 
    
    if (isempty(ylims))
        ylims=[MIN(2) MAX(2)]; 
    end; 

    % Generate means and assign the correct colors
    Data=zeros(size(Y,1),3);
    for i=1:max(c)
        indx=find(c==i);
        Data(:,color(i))=mean(Y(:,indx),2);
    end; 

    % Plot the flatmap 
    caret_plotflatmap_rgb('data',Data,'coord',coord,'topo',topo,...
                    'xlims',xlims,'ylims',ylims,'underlay',underlay,...
                    'dscale',dscale,'threshold',threshold,'border',border); 
    hold on; 
    
    % Superimpose the COGs 
    for q=1:Q
        plot(COG(q,1),COG(q,2),[colorC{c(q)} 'o'],...
            'MarkerSize',markersize,'MarkerFaceColor',colorC{c(q)},...
            'MarkerEdgeColor','k'); 
    end;
    
    if (plotmean) 
        for i=1:max(c)
            indx=find(c==i);
            plot(mean(COG(indx,1)),mean(COG(indx,2)),[colorC{i} 'o'],...
            'MarkerSize',markersize*1.5,'MarkerFaceColor',colorC{i},...
            'MarkerEdgeColor','w');
        end; 
    end; 
    hold off; 
end; 

% Analysis 
switch (test)
    case ''        % No test 
    case 'paired'  % Paired test - assuming same sequence of SN 
            i1=find(c==1); 
            i2=find(c==2);
            if (length(i1)~=length(i2))
                error('for paired test the two classes must have some number of observations.'); 
            end; 
            X=[c' [1:Q/2 1:Q/2]']; 
            hotT2(COG(:,1:2),X,'paired');
            
    case 'independent' 
            X=[c' [1:Q]']; 
            hotT2(COG(:,1:2),X,'independent');
end;         