function S=freesurfer_registerXhem(subj,subject_dir,varargin);

% INPUT:
%   subj: subject name
%   subjects_dir: freesurfer's SUBJECT_DIR
% VARARGIN:
%   'hemisphere',[1 2]  : left / right or both hemispheres
%   'surf_files',{'',''}: Surface files to be resampled
%                       {'.white','.pial','.inflated','.sphere.reg','.sphere'};
%   'curv_files',{'',''}: Curvature files to be resampled
%                       {'.curv','.sulc','.area'};
%   'smoothing',1:        Smoothing iterations applied (otherwise nearest neighbour)
% ---------------------------
% v1.0 Joern Diedrichsen (j.diedrichsen@ucl.ac.uk

%if ~isempty(subject_dir)
%    old_dir=getenv('SUBJECTS_DIR');
%    setenv('SUBJECTS_DIR',subject_dir);
%else
%    subject_dir=getenv('SUBJECTS_DIR');
%end;

subject_dir=getenv('SUBJECTS_DIR');

hemisphere=[1:2]; % Do both hemispheres
target=1;         % Target direction

vararginoptions(varargin,{'hemisphere','target'});

for s=1:length(subj)
    
    subj_dir=[subject_dir filesep subj{s} filesep 'surf'];
    
    if (any(hemisphere==1))
        % Alon: register the subject's left hemisphere to an atlas (fsaverage_sym).        
         % Creates lh.fsaverage_sym.sphere.reg in subject/surf dir. 
        system(['surfreg --s ' subj{s} ' --t fsaverage_sym --lh']);
    end;
    
    if (any(hemisphere==2))
        % Alon: invert right hemisphere
        system(['xhemireg --s ' subj{s}]);
        % Alon: register the subject's inverted right hemisphere to atlas (fsaverage_sym).
        % Creates lh.fsaverage_sym.sphere.reg in subject/xhemi/surf dir.                
        system(['surfreg --s ' subj{s} ' --t fsaverage_sym --lh --xhemi']);
    end;   
end;


%    
%     if (any(hemisphere==1))
%         surf_name=fullfile(subj_dir,'lh.sphere');
%         target=fullfile(subject_dir,'fsaverage_sym','lh.reg.template.tif');
%         outname=fullfile(subj_dir,'lh.fsaverage_sym.sphere.reg');
%         system(['mris_register -curv -annot aparc.annot ' surf_name ' ' target ' ' outname]);
%     end;
%     
%     if (any(hemisphere==2))
%         surf_name=fullfile(subj_dir,'rh.sphere');
%         [V,F]=read_surf(surf_name);
%         V(:,1)=-V(:,1);
%         surf_name=fullfile(subj_dir,'rh.sphere.flipped');
%         freesurfer_write_surf(surf_name,V,F);
%         target=fullfile(subject_dir,'fsaverage_sym','lh.reg.template.tif');
%         outname=fullfile(subj_dir,'lh.fsaverage_sym.X2.sphere.reg');
%         system(['mris_register -curv ' surf_name ' ' target ' ' outname]);
%     end;
%     