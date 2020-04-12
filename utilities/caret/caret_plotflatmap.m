function varargout=caret_plotflatmap(varargin)
% M=caret_plotflatmap
% function M=caret_plotflatmap(varargin)
% Plots vertices as a colormap 
% VARARGIN: 
%   'coord', coord-file structure (flatmap)
%   'topo', topology-file structure (cut)
%   'xlims', Limits of the x-coordinate
%   'ylims', Limits of the y-coordinate
%   'data', Color data to plot (metricfile name, metric structure, or data)
% The colors are not hard-coded, but rely on the current colormap 
% To change the color scheme, just use the command colormap(MAP) 
% For full color control, see caret_plotflatmap_rgb
% j.diedrichsen@ucl.ac.uk, 2011 


border=[]; 
bordersize=5; 
col=1;
M=[];
T=[];
C=[];
topo=[]; 
coord=[];
cscale=[];
interpol='mean'; 
vararginoptions(varargin,{'coord','topo','xlims','ylims','data','col','M','cscale','border','bordersize','interpol'});
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

if (isnumeric(data) | islogical(data))
    D.data=data;
    col=1; 
elseif (isstruct(data))
    D=data;
else 
    D=caret_load(data);
end; 

if (isempty(T))
    T=caret_load(topo);
end;
switch(interpol) 
    case 'mean' 
        d=[D.data(T.data(M.k(:),1),col)+D.data(T.data(M.k(:),2),col)+D.data(T.data(M.k(:),3),col)]'/3;
    case 'mode' 
        d=[D.data(T.data(M.k(:),1),col) D.data(T.data(M.k(:),2),col) D.data(T.data(M.k(:),3),col)]';
        d=mode(d); 
    case 'vert'     % Plot each vertex in it's own color
        d=[D.data(T.data(M.k(:),1),col) D.data(T.data(M.k(:),2),col) D.data(T.data(M.k(:),3),col)]';
end; 
d(isnan(d))=0;      % Kill NaN's 
if (isempty(cscale))
    cscale=[min(d) max(d)];
end;
p=patch(M.X,M.Y,d);

caxis(cscale);
set(gca,'XLim',M.xlims,'YLim',M.ylims);

if (~isempty(border))
    if ~isstruct(border)
        
    end; 
    xB=border.data(:,1); 
    yB=border.data(:,2); 
    i=find(xB>M.xlims(1) & xB<M.xlims(2) & yB>M.ylims(1) & yB<M.ylims(2)); 
    hold on; 
    b=plot(xB(i),yB(i),'k.'); 
    set(b,'MarkerSize',bordersize);
    hold off;
end; 
set(p,'LineStyle','none');
set(p,'EdgeAlpha',0); 

varargout={M,d,p};