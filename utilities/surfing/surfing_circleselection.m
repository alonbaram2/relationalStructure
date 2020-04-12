function [R]=surfing_circleselection(c1,c2,f,circledef,voldef,nodeidxs,linedef,distancemetric)
%   C1,C2:      3xP coordinates for the P vertices: C1 is usually the white, C2 the pial surface  
%   F:          3xQ vertex indices for the Q faces (1-based indices)
%   VOLDEF:     Struct with fields .mat (4x4 voxel to coordinate transformation matrix 
%               from the images to be sampled (1-based)) and .dim (1x3
%               volume dimensions in voxels). Optionally a field .mask with
%               a voxel mask 
%   NODEIDXS:   Sx1 vector for the center vertices for ROIS ([] for all, meaning S==P)
%   LINEDEF:    definition for lines spanning from surface 1 to surface 2 (see SURFING_COORDS2VOXELIDXS)
%   DIST:       Distance metric 'geodesic' (default) or 'euclidian'
%   STEP:       Progress reporting every STEP center vertices (0 for mute) 


% Find coordinates of points on or between the two surfaces.
% Lines between the surfaces are constructed according to LINEDEF.
% ALLCOORDS(I,J,K) is the I-th spatial coordinate (1,2,3)
% for the J-step along the line between the K-th node on surface 1 and 2

allcoords=surfing_nodeidxs2coords(c1, c2,[],linedef);
% Find the voxel indices corresponding to the coordinates of points above
% ALLLINVOXIDXS(I,K) contains the linear index for the voxel
% associated with node I for the K-th step. 
alllinvoxidxs=surfing_coords2linvoxelidxs(allcoords,voldef);
% For each row seperately, duplicates are replaced by zeros.
unqlinvoxidxs=surfing_uniqueidxsperrow(alllinvoxidxs);
clear alllinvoxidxs;
% construct intermediate surface that is the average of the two surfaces
intermediatecoords=squeeze(surfing_nodeidxs2coords(c1, c2,[],[1,0.5,0.5]));
% find the mapping from nodeidxs to the faces that contain the nodes
% this increases the speed of SURFING_SUBSURFACE dramatically
n2f=surfing_nodeidxs2faceidxs(f);

for i= 1:numel(nodeidxs)
    [coordidxs,dist]=surfing_circleROI(intermediatecoords,f,nodeidxs,circledef(1),distancemetric,n2f);
    R{i}= unique(coordidxs(coordidxs>0)); 
    %we could also give back a constant voxel number of the regions with an sec
    %entry in circledef
    %[s, idx]= sort(dist);
    %r=unqlinvoxidxs(coordidxs_S1(idx_S1,:),:);
    %r=unique(r(r>0));
    %R{i}= r(1:circledef(2));
end