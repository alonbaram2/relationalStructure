function S=caret_calcneighbor(S)
% function S=caret_calcneighbor(S)
% finds the neighborhood of each Node, this is a necessary
% computation to find clusters fast an efficiently on the statitical map
% The result is stored in S.Nodes.Neigh, a cell array with numk_nodes
% entries. Each Entry holds the neighbors for that Node
% Remark: notice that the Node-numbering starts at 1, not at 0 as in caret
%---------------------------------------------------
% v.1.1 Joern Diedrichsen 12/04
% jdiedric@bme.jhu.edu
if (isfield(S,'Edges'))
    for n=1:S.num_nodes;
        s1=(S.Edges.data(:,1)==n);
        s2=(S.Edges.data(:,2)==n);
        S.Nodes.Neighbor{n,1}=[S.Edges.data(s1,2);S.Edges.data(s2,1)]';
        if (mod(n,1000)==0)
            fprintf('.',n);
        end;
    end;
    
else
    for n=1:S.num_nodes;
        s=(S.Tiles.data(:,1)==n | S.Tiles.data(:,2)==n | S.Tiles.data(:,3)==n);
        neigh=S.Tiles.data(s,:);
        neigh=neigh(neigh(:)~=n);
        neigh=unique(neigh);
        S.Nodes.Neighbor{n,1}=neigh';
        if (mod(n,1000)==0)
            fprintf('.',n);
        end;
    end;
end;
fprintf('\n',n);
