function Nodes=caret_cluster_pick(S,i,seedvertex);
% function Nodes=caret_cluster_pick(S,i,seedvertex);
% Finds clusters of activations on a surface 
% INPUT: 
%   S: Surface structure from S=caret_getSurface 
%   i: Logical index of the areas included in the set 
%   seedvertex: Vertex which lies in the cluster that you want to extrac 
% OUTPUT: 
%   Nodes: Index of the vertices that you want
% v.2.0 Much more efficient than the old caret_clusters
indx=caret_clusters(S,i); 
if (indx(seedvertex)==0) 
    Nodes=[];       % The seedvertex was outside the excursion set, so return empty 
else 
    Nodes=find(indx==indx(seedvertex)); 
end; 