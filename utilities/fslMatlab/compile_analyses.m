
% compile analyses

clear all

addpath('/vols/Scratch/tmuller/cheese/scripts/NIfTI_20140122')

load('MVPA_svm_1.mat')
AccImage(isnan(AccImage)) = 0;
total_Acc_Image = AccImage;

% total_p = NaN(8, 40, 53373);
% p(isnan(p)) = 0;
% total_p = total_p + p;


for i=2:46
    
    i
    load(['MVPA_svm_' int2str(i) '.mat'])
    
    AccImage(isnan(AccImage)) = 0;
    total_Acc_Image = total_Acc_Image + AccImage;
    
%     p(isnan(p)) = 0;
%     total_p = total_p + p;

end 


save('total_Acc_Image','total_Acc_Image')

nii = make_nii(total_Acc_Image);
save_nii(nii, 'whole_brain.nii.gz');


% 
% load('MVPA_svm_1.mat')
% AccImage(isnan(AccImage)) = 0;
% total_Acc_Image = AccImage;
% 
% load('MVPA_svm_2.mat')
% AccImage(isnan(AccImage)) = 0;
% total_Acc_Image = total_Acc_Image + AccImage;
% 
% load('MVPA_svm_3.mat')
% AccImage(isnan(AccImage)) = 0;
% total_Acc_Image = total_Acc_Image + AccImage;
% 
% load('MVPA_svm_4.mat')
% AccImage(isnan(AccImage)) = 0;
% total_Acc_Image = total_Acc_Image + AccImage;
% 
% load('MVPA_svm_5.mat')
% AccImage(isnan(AccImage)) = 0;
% total_Acc_Image = total_Acc_Image + AccImage;
% 
% load('MVPA_svm_6.mat')
% AccImage(isnan(AccImage)) = 0;
% total_Acc_Image = total_Acc_Image + AccImage;
% 
% load('MVPA_svm_7.mat')
% AccImage(isnan(AccImage)) = 0;
% total_Acc_Image = total_Acc_Image + AccImage;
% 
% load('MVPA_svm_8.mat')
% AccImage(isnan(AccImage)) = 0;
% total_Acc_Image = total_Acc_Image + AccImage;
% 
% load('MVPA_svm_9.mat')
% AccImage(isnan(AccImage)) = 0;
% total_Acc_Image = total_Acc_Image + AccImage;
% 
% load('MVPA_svm_10.mat')
% AccImage(isnan(AccImage)) = 0;
% total_Acc_Image = total_Acc_Image + AccImage;
% 
% load('MVPA_svm_11.mat')
% AccImage(isnan(AccImage)) = 0;
% total_Acc_Image = total_Acc_Image + AccImage;