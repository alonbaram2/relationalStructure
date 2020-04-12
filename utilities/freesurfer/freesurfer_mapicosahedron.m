function S=freesurfer_mapicosahedron(subj,subject_dir,varargin);
% function S=freesurfer_mapicosahedron(subj,subject_dir,varargin);
% Resampels a registered subject surface to a regular isocahedron
% This allows things to happen exactly in atlas space - each vertex number
% corresponds exactly to a anatomical location 
% Makes a new folder, called ['i' subj] that contains the remapped subject
% Uses function mri_surf2surf 
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
% 
curent_dir=pwd;

direct=[getenv('FREESURFER_HOME') filesep 'average' filesep 'surf' ];
old_dir=getenv('SUBJECTS_DIR'); 
setenv('SUBJECTS_DIR',subject_dir); 

structname={'left','right'};
dirname={'LeftHem','RightHem'};
hem={'lh','rh'};
surf_files={'.white','.pial','.inflated','.sphere.reg','.sphere'};
curv_files={'.curv','.sulc','.area'}; 
smoothing=1; 
hemisphere=[1:2]; % Do both hemispheres 

vararginoptions(varargin,{'smoothing','surf_files','curv_files','hemisphere'});

% read freesurfer version
ver_file = fullfile(getenv('FREESURFER_HOME'),'build-stamp.txt');
if exist(ver_file)
    fid     = fopen(ver_file);
    verStr  = fgetl(fid);
    fclose(fid);
    verStr  = strrep(verStr,'-v',' ');
    verStr  = textscan(verStr,'%s%f');
    fsl_ver = verStr{2};
end


% Original subjects dir
orig_dir =[subject_dir filesep subj filesep 'surf'];
new_dir = [subject_dir filesep 'i' subj filesep 'surf'];
if (~exist(new_dir))
    mkdir(new_dir);
end;
cd (orig_dir); 
 
num_surf=length(surf_files); 
files={surf_files{:},curv_files{:}}; 

for h=hemisphere
    [Vico,F]=read_surf([direct filesep 'lh.sphere.reg']);
    
    % Topology file
    for i=1:length(files)
        if i<=num_surf
            [V]=read_surf([orig_dir filesep hem{h} files{i}]);
            
        else
            [V]=read_curv([orig_dir filesep hem{h} files{i}]);
        end;
        B=[]; 
        
        for j=1:size(V,2)
        
            % This is a hack: We transfer the surface files by saving each
            % of the columns (x,y,z) as an individual curvature file. Then
            % transfer them and recompile them together. 
            % Problems may occurr dependign on if mri_surf2surf appends the
            % lh. or rh. in front of temp.curv or not.  TO BE SOLVED. 
            write_curv([orig_dir filesep hem{h} '.temp.curv'],V(:,j),10000);
            system(['mri_surf2surf --srcsubject ' subj ' --hemi ' hem{h} ...
                ' --trgsubject ico  --trgicoorder 7'...
                ' --srcsurfval temp.curv --src_type curv' ...
                ' --trgsurfval tempN.curv --trg_type curv --mapmethod nnf --nsmooth-out ' num2str(smoothing)]);
            %  B(:,i)=read_curv([orig_dir filesep hem{h} '.tempN.curv']);
            % For different alignment, I would have to use -surf_reg option
            % here 
            if fsl_ver <= 5.1
                B(:,j)=read_curv(fullfile(orig_dir,'tempN.curv'));
            else
                % B(:,j)=read_curv(fullfile(orig_dir,'lh.tempN.curv'])); 
                B(:,j)=read_curv(fullfile(subject_dir,subj,'lh.tempN.curv')); 
            end;

%             B(:,j)=read_curv(['tempN.curv']);
        end;
        if i<=num_surf
            freesurfer_write_surf([new_dir filesep hem{h} files{i}],B(:,1:3),F);
        else 
            write_curv([new_dir filesep hem{h} files{i}],B(:,1),size(F,1));
        end; 
    end;
end;
delete temp.curv tempN.curv;
setenv('SUBJECTS_DIR',old_dir); 
cd(current_dir); 