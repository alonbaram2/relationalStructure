function S=freesurfer_mat2registerdat(subj,mean_func,varargin);
% function S=freesurfer_mat2registerdat(subj,mean_func,varargin);
% Makes a register.dat file for a mat to mat alignment from SPM 
% 
% In SPM (note 1-based voxels in SPM, 0 based voxels here!) 
% F0_vox ---> SPACE <--- Ana0_vox
%        M_f        M_a              
% In Freesurfer (0-based voxels) 
% F0_vox ---> FRAS <---------- ARAS <--- Ana0_vox 
%         Tf         register         Ta
% 
% 
% F0_vox = inv(Tf)*R*Ta *A0_vox
% F0_vox =  inv(M_f) * M_a  A0_vox
% Solve for R 
% R = Tf * inv(M_f) * M_a * inv(Ta) 
% see also http://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems?action=AttachFile&do=get&target=fscoordinates.pdf    
% ---------------------------
% v1.0 Joern Diedrichsen (j.diedrichsen@ucl.ac.uk
% 
subj_dir=getenv('SUBJECTS_DIR'); 

Mvox2surf=[-1 0 0 128.5;0 0 1 -128.5;0 -1 0 128.5;0 0 0 1];
anafile=[subj_dir filesep subj filesep 'mri/orig.mgz'];


[status,result] = system(['mri_info ' anafile ' --vox2ras-tkr']); 
A=sscanf(result,'%f'); 
Ta=reshape(A,4,4)';

[status,result] = system(['mri_info ' anafile ' --vox2ras']); 
A=sscanf(result,'%f'); 
Ma=reshape(A,4,4)';

[status,result] = system(['mri_info ' mean_func ' --vox2ras-tkr']); 
A=sscanf(result,'%f'); 
Tf=reshape(A,4,4)';

[status,result] = system(['mri_info ' mean_func ' --vox2ras']); 
A=sscanf(result,'%f'); 
Mf=reshape(A,4,4)';

[status,cres] = system(['mri_info ' mean_func ' --cres']); 
[status,rres] = system(['mri_info ' mean_func ' --rres']); 


R=Tf* inv(Mf) *Ma *inv(Ta);  


fid=fopen('register.dat','w'); 
if (fid==-1) 
    error('could not open register.dat'); 
end; 
fprintf(fid,'%s\n',subj);
fprintf(fid,'%s\n%s\n%2.2f\n',cres,rres,0.15);
fprintf(fid,'%f %f %f\n',R(1,1),R(1,2),R(1,3),R(1,4));
fprintf(fid,'%f %f %f\n',R(2,1),R(2,2),R(2,3),R(2,4));
fprintf(fid,'%f %f %f\n',R(3,1),R(3,2),R(3,3),R(3,4));
fprintf(fid,'%f %f %f\n',R(4,1),R(4,2),R(4,3),R(4,4));
fprintf(fid,'round\n'); 
fclose(fid);