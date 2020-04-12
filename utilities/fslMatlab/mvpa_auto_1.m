% roving searchlight MVPA w/ OneVsAll SVM

clear ; close all; clc
run=1;

% mask=read_avw('gm_mask.nii.gz');
% mask = padarray(mask,[5,5,5]);
% mask = load('mask.mat');
% mask = mask.mask;
mask=padarray(read_avw(['mask.nii.gz']),[5,5,5]); 
ix=find(mask>0);
[x,y,z]=ind2sub(size(mask),ix);
%nSearchlights=length(x);
nSearchlights=5000;
s1 = (run-1)*nSearchlights + 1;
s2 = min((run)*nSearchlights, length(ix));
nSearchlights=s2-s1 + 1; % prevents breaking

nExemplars=40;
nClasses = 8;

fprintf(['Running MVPA for ' int2str(nSearchlights) ' searchlights!\n']);

%nSearchlights=10;

r = 3; % radius (in voxels) of each searchlight

[xcoords,ycoords,zcoords]=ndgrid(1:size(mask,1), 1:size(mask,2), 1:size(mask,3));

% generate matrix with indices of searchlights
fprintf(['Finding searchlight indices\n']);
%for i=1:nSearchlights
for i=s1:s2
    j=(i-s1)+1;
    sphere_ix(:,j) = find(sqrt(((xcoords-x(i)).^2) + ((ycoords-y(i)).^2) + ((zcoords-z(i)).^2))<r);
end


% make a huge matrix of: examples x voxels (=features) x searchlights
fprintf(['Reading in data\n']);
X=NaN(nExemplars, size(sphere_ix,1), nSearchlights);
for m = 1:nExemplars; % no. of examples
    fprintf(['Reading COPE ' int2str(m) '\n']);
    temp=padarray(read_avw(['pe' int2str(m) '.nii.gz']),[5,5,5]);
    for i = 1:nSearchlights  % no. of searchlights/loop over every voxel
        X(m,:,i)=temp(sphere_ix(:,i))';
    end;
end;

X = X-repmat(mean(X,1),40,1); X = X./repmat(std(X),40,1);
X_not_nan = squeeze(sum(~isnan(X(1,:,:)),2));
X(isnan(X))=0;
classes = dlmread('y');
% 
  % now run MVPA on each of the layers of matrix X (i.e. each searchlight)
fprintf(['Running SVM\n']);
S=NaN(nExemplars, nClasses, nSearchlights);
L=NaN(nExemplars, nClasses, nSearchlights);
for i=1:nSearchlights %size(X,3)
%       j=(i-s1)+1;
      fprintf(['Searchlight ' int2str(i) '\n']);
      [S(:,:,i), L(:,:,i)] = mvpa_svm_OneVsAll(X(:,:,i), classes);
      
      if mod(i,100)==0 
        dlmwrite('S',S);
      end
end

save([ 'MVPA_svm_' int2str(run)]); 

p=1./(1+exp(-S));
p=permute(p,[2 1 3]);
for i=1:nSearchlights
   p(:,:,i)=p(:,:,i)./repmat(sum(p(:,:,i)),8,1); % make probabilities add up to 1 across classes for each exemplar
end
%plot_p;
for i=1:40
    pShift(:,i,:) = circshift(p(:,i,:),4-classes(i),1);
end

AccVals = squeeze(mean(pShift(4,:,:),2));
AccImage = zeros(size(mask));
AccImage(ix(s1:s2)) = AccVals;
AccImage = AccImage(6:end-5, 6:end-5, 6:end-5); % unpad image to get back to MNI space
save(['AccImage' int2str(run) '.txt'],'AccImage');
%save_avw(AccImage,'Acc.nii.gz','f',[3 3 3 3])

save([ 'MVPA_svm_' int2str(run)]);