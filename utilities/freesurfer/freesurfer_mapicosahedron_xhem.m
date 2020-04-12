function S=freesurfer_mapicosahedron_xhem(subj,subject_dir,varargin);
% function S=freesurfer_mapicosahedron(subj,subject_dir,varargin);
% Resampels a registered subject surface to a regular isocahedron
% This allows things to happen exactly in atlas space - each vertex number
% corresponds exactly to a anatomical location 
% Makes a new folder, called ['x' subj] that contains the remapped subject
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
% 5/Sep/13 Naveed: added a check for the freesurfer version to determine which 
%                  of the tempN.curv or lh.tempN.curv is generated for xhem registration
current_dir=pwd;

direct=[getenv('FREESURFER_HOME') filesep 'average' filesep 'surf' ];
old_dir=getenv('SUBJECTS_DIR'); 
% if (isempty(subject_dir))
%    subject_dir=getenv('SUBJECTS_DIR'); 
% else 
%    setenv('SUBJECTS_DIR',subject_dir); 
% end;
subject_dir=getenv('SUBJECTS_DIR');
structname={'left','right'};
dirname={'LeftHem','RightHem'};
hem={'lh','rh'};
surf_files={'.white','.pial','.inflated'};
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

    
new_dir = [subject_dir filesep 'x' subj filesep 'surf'];
if (~exist(new_dir))
    mkdir(new_dir);
end;

num_surf=length(surf_files); 
files={surf_files{:},curv_files{:}}; 

for h=hemisphere
    % Original subjects dir
    if (h==1) 
        orig_dir =fullfile(subject_dir,subj,'surf');
    else
        orig_dir =fullfile(subject_dir,subj,'xhemi','surf');
    end; 
    cd(orig_dir); 
    % [vertex_coords, faces]
    [Vico,F]=read_surf([direct filesep 'lh.sphere.reg']); % /opt/fmrib/FreeSurfer_releases/6.0/average/surf/lh.sphere.reg

    
    % Topology file
    for i=1:length(files)
        if i<=num_surf
            % e.g /home/fs0/abaram/scratch/freesurfer-subjects/s99/surf/lh.white
            [V]=read_surf([orig_dir filesep 'lh' files{i}]);
            
        else
            [V]=read_curv([orig_dir filesep 'lh' files{i}]);
        end;
        B=[]; 
        for j=1:size(V,2) % Alon: number of coords - usually 3
            % Alon: write the coords of the vertices (one dimension at a
            % time) as curv files (associate verteces with a single numer)
            write_curv([orig_dir filesep 'lh.temp.curv'],V(:,j),10000); % e.g [$SUBJECTS_DIR/s99/surf/lh.temp.curv, only x coords of V, number of faces (why only 10000?! I think it's just not used)]
            if (h==1)
                system(['DYLD_LIBRARY_PATH=/opt/fmrib/FreeSurfer_releases/6.0/lib/gcc/lib '...
                    'mri_surf2surf --srcsubject ' subj ' --hemi lh' ...
                    ' --trgsubject ico  --trgicoorder 7'...
                    ' --surfreg fsaverage_sym.sphere.reg'...
                    ' --srcsurfval temp.curv --src_type curv' ...
                    ' --trgsurfval tempN.curv --trg_type curv --mapmethod nnf --nsmooth-out ' num2str(smoothing)]);
            else 
                % Reslice the left hemisphere of the xhemi subject - this
                % is originally the right hemisphere
                system(['DYLD_LIBRARY_PATH=/opt/fmrib/FreeSurfer_releases/6.0/lib/gcc/lib '...
                    'mri_surf2surf --srcsubject ' subj '/xhemi --hemi lh' ...
                    ' --trgsubject ico  --trgicoorder 7'...
                    ' --surfreg fsaverage_sym.sphere.reg'...
                    ' --srcsurfval temp.curv --src_type curv' ...
                    ' --trgsurfval tempN.curv --trg_type curv --mapmethod nnf --nsmooth-out ' num2str(smoothing)]);
            end; 
%              B(:,i)=read_curv([orig_dir filesep hem{h} '.tempN.curv']);
            %  For different alignment, I would have to use -surf_reg option
            %  here 
            
            % Alon: read in the resampled curv file of a single dim of the
            % coordinates per vertex. At the end of the j loop, B will have
            % the all 3 resamped dims of the coordinates. 
            if fsl_ver <= 5.1
                B(:,j)=read_curv(fullfile(orig_dir,'tempN.curv'));
            else
                B(:,j)=read_curv(fullfile([orig_dir filesep 'lh.tempN.curv']));
            end;
        end;
        % write the surfaces (3 col) and curvs (1 col).
        if i<=num_surf
            freesurfer_write_surf([new_dir filesep hem{h} files{i}],B(:,1:3),F);
        else 
            write_curv([new_dir filesep hem{h} files{i}],B(:,1),size(F,1));
        end; 
    end;
    
    % deleting temporary files
    try
        delete(fullfile(orig_dir,['lh.temp.curv']));
        if fsl_ver <= 5.1 
            delete(fullfile(orig_dir,['tempN.curv']));
        else
            delete(fullfile(orig_dir,['lh.tempN.curv']));
        end
    catch err
        disp(getReport(err));
    end;
    
end;
setenv('SUBJECTS_DIR',old_dir); 