function [A,s]=caret_calc_EC(S,mask)
% function [A,stats]=caret_calc_EC(S,mask)
% Calculates the area, circumference and Euler Characteristic of a subset of a surface
% relies on area, edge-length from caret_calcarea
% INPUT :
%   S: surface (caret_getsurface) 
%   mask: logical vector of length of Nodes 
% OUTPUT: 
%   A=[Euler characteristic circference surface]; 
%   s: Stucture with basic information on search region 
%---------------------------------------------------
% v 2.0 Joern Diedrichsen 12/20
% v.2.1 Small fixes on 
% j.diedrichsen@ucl.ac.uk

% --------------------------------------------------
% First find the right Tiles
Tmask=zeros(S.num_tiles,3);
for j=1:3
    Tmask(:,j)=mask(S.Tiles.data(:,j));
end;
Tmask=all(Tmask,2);

% --------------------------------------------------
% Then recalculate the edges:
D=S.Tiles.data(Tmask,:);
L=S.Tiles.edgelength(Tmask,:);

Edges.data=[];
for i=1:3
    i1=D(:,i);i2=D(:,mod(i,3)+1);
    Edges.data=[Edges.data;[i1 i2]];
end;
Edges.data=sort(Edges.data,2);
[Edges.data,Edges.indx,point]=unique(Edges.data,'rows');
Edges.num_occ=zeros(length(Edges.data),1);
for i=1:length(point)
    Edges.num_occ(point(i))=Edges.num_occ(point(i))+1;
end;
Edges.Length=L(Edges.indx); % Copy edges from tiles 

s.num_nodes=sum(mask>0);
s.num_tiles=sum(Tmask); 
s.num_edges=size(Edges.data,1);
s.bad_tiles=sum(S.Tiles.area==0);
s.bad_nodes=sum(S.Nodes.area==0);
s.bad_edges=sum(S.Edges.Length==0);
s.triple_edges=sum(S.Edges.num_occ>2); 
s.single_edges=sum(S.Edges.num_occ==1);

%------------------------------------------------------
% Now compute the EC, half the perimeter and surface area of the surface
% Make sure that any nodes that are totally disconnected are not counted
% also make sure that there are no bad tiles
A(1)=sum(mask)-s.num_edges+sum(Tmask);
A(2)=sum(Edges.Length(Edges.num_occ==1))/2;
A(3)=sum(S.Nodes.area(logical(mask)));