function varargout=caret_plotflatmap_rgb(varargin)
% varargout=caret_plotflatmap_rgb(varargin)
% Plots a surface patch as a flatmap - using a functional RGB overlay
% VARARGIN:
%   'coord'     : Name of coordinate file
%   'topo'      : Name of topo file
%   'data'      : Data to plot. Its either a (Vx3) value, which uses a
%   'xlims'     :
%   'ylims'     :
%   'underlay'  :
%   'uscale'    : [min max] Mapping black to white for the underlay
%   'M'         : Surface structure for faster plotting
%   'dscale'    : [min max] Mapping black to full rgb strength
%   'threshold' : Display threshold for overlay (by default -inf)
%   'border'    : Border file or structure to plot
%   'alpha'     : Alpha value for overlay
% OUTPUT:
%    M          : Summary structure of surface patch (returned for speedier
%                   subsequent call
%    d          :
%    p          :
border=[];
col=1;
M=[];
T=[];
C=[];
topo=[];
coord=[];
cscale=[];
underlay=[];
uscale=[];
umap = gray;
data=[];
dscale=[];
dmap = jet;
alpha=1;
threshold = -inf;

vararginoptions(varargin,{'coord','topo','xlims','ylims','data',...
    'underlay','uscale','umap','M',...
    'dscale','dmap','border','alpha','threshold'});

if (isempty(xlims))
    error('Must supply xlims');
end;
if (isempty(ylims))
    error('Must supply xlims');
end;

if (isempty(M))
    C=caret_load(coord);
    T=caret_load(topo);
    
    % Get the X-Y coordinates for all tiles
    X(1,:)=C.data(T.data(:,1),1);
    X(2,:)=C.data(T.data(:,2),1);
    X(3,:)=C.data(T.data(:,3),1);
    Y(1,:)=C.data(T.data(:,1),2);
    Y(2,:)=C.data(T.data(:,2),2);
    Y(3,:)=C.data(T.data(:,3),2);
    
    % Find all tiles that have a single vertex (or more) in the image
    M.k=find(any(X>xlims(1) & X<xlims(2),1) & any(Y>ylims(1) & Y<ylims(2),1));
    M.X=X(:,M.k);
    M.Y=Y(:,M.k);
    
    M.xlims=xlims;
    M.ylims=ylims;
end;

if (isempty(T))
    T=caret_load(topo);
end;

% Generate mapping for the underlay
% Determine the shading of the faces by the vertices and scale the color assignment
if (size(underlay,2)==1)                            % Single value
    u=[underlay(T.data(M.k(:),1),1) underlay(T.data(M.k(:),2),1) underlay(T.data(M.k(:),3),1)];
    numcols=size(umap,1); % Size of underlay colormap
    if (isempty(uscale))
        uscale=[min(underlay) max(underlay)];
    end;
    uindx=ceil((u-uscale(1))./(uscale(2)-uscale(1))*numcols);
    uindx(uindx<1)=1;
    uindx(uindx>numcols)=numcols;
    
    % Now assign the RGB color to each face
    for i=1:3
        for j=1:size(uindx,1)
            uCOL(j,:,i) = umap(uindx(j,:),i);
        end;
    end;
elseif (size(underlay,2)==3)                       % Already RGB values: must scale
    for i=1:3
        u(:,:,i)=[underlay(T.data(M.k(:),1),i) underlay(T.data(M.k(:),2),i) underlay(T.data(M.k(:),3),i)];
    end;
    if (isempty(uscale))
        uscale=[min(u) max(u)];
    end;
    uCOL=(u-uscale(1))./(uscale(2)-uscale(1));
    uCOL(uCOL<0)=0;
    uCOL(uCOL>1)=1;
else                                                % Empty:
    uCOL = ones(length(M.k(:)),3,3);
end;

% Generate the mapping for the overlay
if (size(data,2)==1)      % Single value
    d=[data(T.data(M.k(:),1),1) data(T.data(M.k(:),2),1) data(T.data(M.k(:),3),1)];
    numcols=size(dmap,1); % Size of underlay colormap
    if (isempty(dscale))
        dscale=[min(d(:)) max(d(:))];
    end;
    dindx=ceil((d-dscale(1))./(dscale(2)-dscale(1))*numcols);
    dindx(dindx<1 | isnan(dindx))=1;  %
    dindx(dindx>numcols)=numcols;
    % Now assign the RGB color to each face
    for i=1:3
        for j=1:size(dindx,1)
            dCOL(j,:,i) = dmap(dindx(j,:),i);
        end;
    end;
    
    % Set value below threshold to NaN
    for i=1:3
        for j=1:3
            dCOL(d(:,j)<threshold | isnan(d(:,j)),j,i)=NaN;
        end;
    end;
    
    
elseif (size(data,2)==3)                        % RGB colors
    for i=1:3
        d(:,:,i)=[data(T.data(M.k(:),1),i) data(T.data(M.k(:),2),i) data(T.data(M.k(:),3),i)];
    end;
    
    
    if (isempty(dscale))
        dscale=[min(d,[],2) max(d,[],2)];
    else
        [i,j]=size(dscale);
        if (i==2 && j==3)
            dscale=dscale';
        elseif (i==3 && j==2)  %Correct
        elseif (i==1 && j==2)
            dscale=[dscale;dscale;dscale];
        else
            error('dscale must be 1x2, 2x3 or 3x2 matrix');
        end;
    end;
    if (length(threshold)==1)
        threshold=ones(1,3)*threshold;
    end;
    
    for i=1:3
        range=dscale(i,2)-dscale(i,1);
        if range<0
            error('dscale max needs to be bigger than dscale min');
        end;
        if (range>0)
            dCOL(:,:,i)=(d(:,:,i)-dscale(i,1))/range;
        else
            dCOL(:,:,i)=d(:,:,i);
        end;
    end;
    
    % Limit to maximum color range 
    dCOL(dCOL<0)=0; 
    dCOL(dCOL>1)=1; 
    
    
    % Set value below threshold to NaN
    for j=1:3   % Nodes of triangle 
        indx=all(bsxfun(@lt,squeeze(d(:,j,:)),threshold),2);  % Where are all under threshold? 
       	dCOL(indx,j,:)=NaN; 
    end; 

else
    dCOL = nan(length(M.k(:)),3,3);
end;

% Find empty places on the RGB and replace with underlay
empty_nodes=(isnan(dCOL));

% Blend the underlay with overlay:
COL=(1-alpha)*uCOL+alpha*dCOL;
COL(empty_nodes)=uCOL(empty_nodes);

% Draw the patches
p=patch(M.X,M.Y,permute(COL,[2,1,3]));
set(gca,'XLim',M.xlims,'YLim',M.ylims);
set(p,'LineStyle','none');
set(p,'EdgeAlpha',0);

% Draw border
if (~isempty(border))
    for i=1:length(border)
        xB=border(i).data(:,1);
        yB=border(i).data(:,2);
        i=find(xB>M.xlims(1) & xB<M.xlims(2) & yB>M.ylims(1) & yB<M.ylims(2));
        hold on;
        plot(xB(i),yB(i),'k.');
    end;
    hold off;
end;

varargout={M,p};