function caret_mirrordata(varargin);
% Mirrors data between the left and the right hemisphere
% Needs the alignment to fsaverage_sym for accurate result yet 

% Test how asymmetric surface shape for the fsaverage compared to
% fsaverage_sym is 
% Note that the spheres are exactly the same for the left and the 
% Right hemispheres. 


% Also, the surface-shape files a
VL=caret_load('fsaverage_sym/LeftHem/lh.SPHERE.REG.coord');
VR=caret_load('fsaverage_sym/RightHem/rh.SPHERE.REG.coord');
VR.data(:,1)=-VR.data(:,1); % Flip the sphere 
match=zeros(VL.num_nodes,1); 
indx=zeros(VL.num_nodes,1); 

% Now seek the corresponding vertices
for i=1:VL.num_nodes
    A=bsxfun(@minus,VR.data,VL.data(i,:)); 
    [match(i,1),indx(i,1)]=min(sum(A.^2,2)); 
    if (mod(i,1000)==0)
        fprintf('.'); 
    end;
end; 

save mirrorSPHERE.mat match indx 
