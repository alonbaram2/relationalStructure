function S=freesurfer_reconall(subject_dir,subj_name,anatomical,varargin);
% function S=freesurfer_reconall
% Simply call recon all
% ---------------------------
% v1.0 Joern Diedrichsen (j.diedrichsen@ucl.ac.uk
%

if ~isempty(subject_dir)
    old_dir=getenv('SUBJECTS_DIR');
    setenv('SUBJECTS_DIR',subject_dir);
else
    subject_dir=getenv('SUBJECTS_DIR');
end;

xhemi=1; 

vararginoptions(varargin,{'hemisphere','xhemi'});

system(['recon-all -s ' subj_name ' -i ' anatomical ' -all -cw256']);

