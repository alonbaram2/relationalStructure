function S=caret_calcarea(S,flag)
% Calculates the area, circumference and Euler Characteristic of an surface
% The function excludes Nodes and Tiles that do not contain any surface
% area.
% The basic values of the area and length of Nodes and Tiles is stored in
% the respective fields for later fast calculation of cluster-sizes and
% corrected p-values
%---------------------------------------------------
% v.1.1 Joern Diedrichsen 12/04
% v 2.0: Joern Diedrichsen 12/10
if(nargin<2)
    flag='';
end;
%---------------------------------------------------
% If necessary calculate area of each triangle and keep track of edges for later
if (~isfield(S.Tiles,'area') | ~isfield(S.Tiles,'edgelength') | strcmp(flag,'force'))
    % fprintf('Calculating surface of Tiles\n');
    for i=1:3
        i1=S.Tiles.data(:,i);i2=S.Tiles.data(:,mod(i,3)+1);
        x=S.Nodes.data(i1,:)-S.Nodes.data(i2,:);
        S.Tiles.edgelength(:,i)=sum(x.^2,2);
    end;
    x=S.Tiles.edgelength;
    S.Tiles.edgelength=sqrt(S.Tiles.edgelength);
    S.Tiles.area=0.25*sqrt(2.*x(:,1).*x(:,2)+2.*x(:,2).*x(:,3)+2.*x(:,1).*x(:,3)-sum(x.^2,2));
    indx=find(imag(S.Tiles.area));
    S.Tiles.area(indx)=0;
end;
%------------------------------------------------------
% Calculate area of each node
if (~isfield(S.Nodes,'area')  | strcmp(flag,'force'))
    % fprintf('Calculating surface of Nodes\n');
    S.Nodes.area=zeros(S.num_nodes,1);
    for i=1:3
        indx=S.Tiles.data(:,i);
        for j=1:length(indx)
            S.Nodes.area(indx(j))=S.Nodes.area(indx(j))+S.Tiles.area(j)/3;
        end;
    end;
end;

%------------------------------------------------------
% Find the unique Edges from tiles
if (~isfield(S, 'Edges')| ~isfield(S, 'Length')  | strcmp(flag,'force'))
    fprintf('Calculating Edges\n');
    S.Edges.data=[];
    for i=1:3
        i1=S.Tiles.data(:,i);i2=S.Tiles.data(:,mod(i,3)+1);
        S.Edges.data=[S.Edges.data;[i1 i2]];
    end;
    S.Edges.data=sort(S.Edges.data')';
    [S.Edges.data,S.Edges.indx,point]=unique(S.Edges.data,'rows');
    S.Edges.num_occ=zeros(length(S.Edges.data),1);
    for i=1:length(point)
        S.Edges.num_occ(point(i))=S.Edges.num_occ(point(i))+1;
    end;
    S.Edges.Length=S.Tiles.edgelength(S.Edges.indx);
    S.num_edges=size(S.Edges.data,1);
end;

%------------------------------------------------------
% Now compute the EC, half the perimeter and surface area of the surface
% Make sure that any nodes that are totally disconnected are not counted
% also make sure that there are no bad tiles
S.bad_tiles=sum(S.Tiles.area==0);
S.bad_nodes=sum(S.Nodes.area==0);
S.bad_edges=sum(S.Edges.Length==0);
S.triple_edges=sum(S.Edges.num_occ>2);
S.single_edges=sum(S.Edges.num_occ==1);
S.A(1)=S.num_nodes-S.num_edges+S.num_tiles;
S.A(2)=sum(S.Edges.Length(S.Edges.num_occ==1))/2;
S.A(3)=sum(S.Nodes.area);