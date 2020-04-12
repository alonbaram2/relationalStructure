function collapseRdmXRun(rootData,sub,glm,sl,nElem)

%%  load the expanded (2*nElem x 2*Elem) RDM
rdmDir = fullfile(rootData,'RDMs',glm,sub);
dataFile = fullfile(rdmDir,['dist_xRun_' sl '.nii']);
outputFile = fullfile(rdmDir,['dist_xRun_' sl '_collapsed.nii']);

% get NIFTI header info
V = spm_vol(dataFile);
% read NIFTI array
dist_nii = read_avw(dataFile);

%%

% find (flattened) indeces of lower left nElem x nElem square of the expanded 
% 2*nElem x 2*nElem RDM. This can be thought of as a non-symetrical nElem x nElem
% RDM.
indMat = logical(squareform((blkdiag(ones(nElem),ones(nElem))-1)*-1));
Xvec = dist_nii(:,:,:,indMat); % the non-symetrical RDM, flattened
% Change to a square matrix, instead of flatened matrix.
Xmat = reshape(Xvec,size(Xvec,1),size(Xvec,2),size(Xvec,3),nElem,nElem);
% Average lower and upper triangles - i.e. average the two distances
% (d1,d2) between each couple of conditions: d1 is between condition i in run 
% 1 and condition j in run 2 and d2 is between conditino i in run 2 and
% condition j in run 1) 
Xmat = (Xmat + permute(Xmat,[1,2,3,5,4])) / 2;
% Extract diagonal elements (distance between a condition to itself - across 
% runs). We will later append these elements to the flattned version of the 
% symmetric RDM of the off-diagonal elements
indDiag = 1:nElem+1:nElem*nElem;
diagDist = Xmat(:,:,:,indDiag);
% zero diagonal (to make it a distance matrix which can be easily
% flattened).
Xmat = reshape(Xmat,size(Xmat,1),size(Xmat,2),size(Xmat,3),nElem*nElem);
Xmat(:,:,:,indDiag) = 0;
Xmat = reshape(Xmat,size(Xmat,1),size(Xmat,2),size(Xmat,3),nElem,nElem);

% For each voxel, flatten the (off-diagonals) RDM using squareform. This is quite
% unefficient - there's probably a vectorised way of doing this. 
xRunsDist = nan(size(Xmat,1),size(Xmat,2),size(Xmat,3),nElem*(nElem-1)/2);
for i=1:size(Xmat,1)
    for j=1:size(Xmat,2)
        for k=1:size(Xmat,3)
            xRunsDist(i,j,k,:) = squareform(reshape((Xmat(i,j,k,:,:)),nElem,nElem));
        end
    end
end
% append the diagonal elemnts at the end of the flattened off-diagonal
% vector. 
nii = cat(4,xRunsDist,diagDist);
% save
save4Dnii(nii,V,outputFile);