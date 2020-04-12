function index=caret_clusters(S,i);
% function index=caret_cluster_dmperm(S,i);
% Finds clusters of activations on a surface 
% INPUT: 
%   S: Surface structure from S=caret_getSurface 
%   i: Logical index of the areas included in the set 
% OUTPUT: 
%   index: Index of cluster number 
% v.2.0 Much more efficient than the old caret_clusters
G=sparse(S.Edges.data(:,1),S.Edges.data(:,2),ones(size(S.Edges.data,1),1),S.num_nodes,S.num_nodes);
index=double(i);
in=find(i); 
G=G(i,i);
G=G+G'+speye(size(G,1)); 
[p,q,r,s]=dmperm(G); 
cl_size=r(2:end)-r(1:end-1); 
for c=1:length(r)-1 
    index(in(p(r(c):r(c+1)-1)))=c; 
end;  
