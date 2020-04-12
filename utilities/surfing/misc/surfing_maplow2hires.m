function [low2hi,ds]=surfing_maplow2hires(clow, chigh, epsilon)
% finds mapping between nodes from low and high res surfaces
%
% [L2H,D]=SURFING_MAPLOW2HIRES(CLOW,CHIGH[,EPSILON]) 
% INPUTS:
%   CLOW:     3xNLOW  coordinates for N nodes from low  resolution surface
%   CHIGH:    3xNHIGH "                          " high "                "
% OPTIONAL INPUT:
%   EPSILON:  Scalar value of maximum difference between nodes in CLOW and
%             CHIGH. If omitted, EPSILON=1e-4 and surfaces are assumed to
%             be generated using mapicosahedron (see below)
% OUTPUTS:
%   L2H:    NLOWx1 node mapping from low to high resolution surface nodes.
%           If L2H(I)==J, then node J on CHIGH is the nearest node to node
%           J on CLOW.
%   D:      NLOWx1 euldiain distance from node I in CLOW and node D(I) in
%           CHIGH.
%
% If EPSILON is omitted, it is assumed that these surfaces are created using 
% AFNI's or Freesurfer's MapIcosehedron, with the number of linear divisions 
% (LD) of the "large" edges of the original edges being LD1=K for CLOW and 
% LD2=P*K for CHIGH (where P is an integer). 
% 
% NNO Aug 2009, updated May 2010.

frommapico=nargin<3 || isempty(epsilon); % assume from mapicosahedron if epsilon is omitted

% number of vertices
if size(clow,1) ~= 3 || size(chigh,1) ~= 3, error('Coordinates should be 3xN'); end
    
nlow=size(clow,2);
nhi=size(chigh,2);

if nlow>=nhi
    warning('Lo res surface does not have less vertices than hi res surface; check your input');
end

if frommapico 
    % do a couple of checks
    plow=sqrt((nlow-2)/10);
    phi=sqrt((nhi-2)/10);

    if (floor(plow) ~= plow || floor(phi) ~= phi)
        warning('Surfaces do not seem to be created with mapicosehedron; will continue anyway');
    end

    k=phi/plow;

    if frommapico && abs(round(k)- k)>.001
        error('Number of subdivisions do not match: hi res should be multiple of low res, but fraction is %d', k);
    end
end

low2hi=zeros(nlow,1);
ds=zeros(nlow,1);
for k=1:nlow
    % distance from low res vertex k to all high res vertices
    hids=surfing_eucldist(clow(:,k),chigh);
    
    % find nearest high res vertex
    [d,idx]=min(hids);
    
    
    if frommapico
        if d>1e-4
            error('Low res node %d has smallest distance %.3f > 0.001 to high res node %d, this is more than expected from mapicosehdron',k,d,idx);
        end
    elseif d>epsilon
        error('Low res node %d has smallest distance %.3f to high res node %d, exceeding epsilon=%.3f',k,d,idx,epsilon);
    end
       
    low2hi(k)=idx;
    ds(k)=d;
end

