function nbrs=surfing_surface_nbrs(f)
% Finds the neighbours on a surface
%
% INPUT:
%  F      Nx3 face matrix for N faces
% OUTPUT:
%  NRBS   NxP neighbour matrix, if each node has at most D neighbours
%         The K-th row contains the neighbours of node K, and values zero
%         if node K has less than P neighbours
%
% NNO Sep 2011

[n,k]=size(f);
if n==3 && k~=3
    warning('surfing:surface_nbrs','surfing_nbrs expected Nx3 face matrix; will transpose ');
    f=f';
    [n,k]=size(f);
end

[fi,nfi]=surfing_invertmapping(f); % vertices to facs
[nv,mx]=size(fi); % number of vertices, and max number of neighbours per face
allnbrs=zeros(nv,mx*k); % vertices contained in the faces containing the center vertex
for j=1:mx
    msk=nfi>=j; % omit zero values
    allnbrs(msk,(j-1)*k+(1:k))=f(fi(msk,j),:);
end

unqnbrs=surfing_uniqueidxsperrow(allnbrs);
nbrs=zeros(size(unqnbrs)-[0 1]);
for j=1:nv
    % remove self and empty values
    msk=unqnbrs(j,:) ~= 0 & unqnbrs(j,:) ~= j;
    nbrs(j,1:sum(msk))=unqnbrs(j,msk);
end






