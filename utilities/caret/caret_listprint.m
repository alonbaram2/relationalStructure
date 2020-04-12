function caret_listprint(style,TabDat);
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
%   Node : Number of Node at Maximum
%-----------------------------------------------------------------------
% v.2.0. Joern Diedrichsen 12/4/10 j.diedrichsen@ucl.ac.uk
switch(style)
    % -----------------------------------------------------------------
    % -Print ASCII text table
    case 'txtlist'
        
        %-Table Title
        %-----------------------------------------------------------------------
        fprintf('\n\nSTATISTICS: %s\n',TabDat.tit);
        fprintf('%c','='*ones(1,80)), fprintf('\n');
        
        %-Table header
        %-----------------------------------------------------------------------
        fprintf('%s\t',TabDat.hdr{1:end-1});fprintf('%s\n',TabDat.hdr{end});
        fprintf('%c','-'*ones(1,80)), fprintf('\n');
        
        
        %-Table data
        %-----------------------------------------------------------------------
        for i = 1:size(TabDat.dat,1)
            for j=1:size(TabDat.dat,2)
                fprintf(TabDat.fmt{j},TabDat.dat{i,j});
                fprintf('\t');
            end
            fprintf('\n');
        end
        %      fprintf('%s\n',TabDat.str)
        fprintf('%c','-'*ones(1,80)), fprintf('\n');
        
        %-Table footer
        %-----------------------------------------------------------------------
        fprintf('%s\n',TabDat.ftr{:});
        fprintf('%c','='*ones(1,80)), fprintf('\n\n');
end;
