function Msurf2space=caret_importfreesurfer(subj,freesurfer_dir,caret_dir,varargin);
% Imports a surface reconstruction from Freesurfer automatically into caret
% Makes a new spec file and moves the coord-files to respond to World,
% rather than to RAS_tkreg coordinates.
% see also http://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems?action=AttachFile&do=get&target=fscoordinates.pdf
% Expects a particular directory structure (see below)
% IMPORT:
%   subj: Name of Subject to be transformed
%   Freesurfer_dir: $SUBJECTS_DIR from freesurfer. Files are in
%                   freesurfer_dir/subj/surf
%   Caret_dir:      Caret Output directory. Files will be put into
%                   caret_dir/subj/LeftHem and caret_dir/subj/RightHem
% VARARGIN:
%   'surfaceStr',{'',''}  :Names of Surfaces to be imported
%   'surftype',{'',''}    : Caret Surface types associated with this
%   'shape',{'',''}       : Names of files to be put into the surfaceshape file
%   'shapesign',[]        : Scale value for each column
% ---------------------------
% v1.0 Joern Diedrichsen (j.diedrichsen@ucl.ac.uk
%

structname={'left','right'};
dirname={'LeftHem','RightHem'};
hem={'lh','rh'};
surfaceStr={'white','pial','inflated','sphere','sphere.reg'};
surftype={'FIDUCIAL','OTHER','INFLATED','SPHERICAL','SPHERICAL'};
shape={'curv','sulc','area'}; % 'area'
shapesign=[-1 -1 1];
noshift=0;

old_dir=getenv('SUBJECTS_DIR');
setenv('SUBJECTS_DIR',freesurfer_dir);

if (subj(1)=='i') % Subject is icosahedron - look up original subject
    subjO=subj(2:end);
    flipR=0;
elseif (subj(1)=='x') % Subjects is cross-hemisphere aligned - flip right hemisphere back
    subjO=subj(2:end);
    flipR=1;
else
    subjO=subj;
    flipR=0;
end;

vararginoptions(varargin,{'hem','surfaceStr','surftype','shape','subjO','flipR'},{'noshift'});

surf_dir=[freesurfer_dir filesep subj filesep 'surf'];
for h=1:length(hem)
    % Make spec file
    SPEC.header{1}='BeginHeader';
    SEPC.header{2}='category Individual';
    SPEC.header{3}='encoding ASCII';
    SPEC.header{4}='space OTHER';
    SPEC.header{5}='species Human';
    SPEC.header{6}=['structure' structname{h}];
    SPEC.header{6}=['subjects' subj];
    SPEC.header{7}=['EndHeader'];
    SPEC.data={};
    
    outd=([caret_dir filesep subj filesep dirname{h}]);
    if (~exist(outd,'dir'))
        mkdir(outd);
    end;
    
    % Figure out the shifting of coordinate systems:
    % Freesurfer uses vertex coordinates in respect to
    % the center of the 256x256x256 image.
    % Independent of the real zero point in the original image
    % So to find a transform of the
    % Mvox2surf: Transform of voxels in 256x256 image to surface vertices
    % Mvox2space: Transform of voxel to subject space
    anafile=[freesurfer_dir filesep subjO filesep 'mri/brain.mgz'];
    
    % Alon: Get transform conformed2subjectSurface
    [status,result] = system(['DYLD_LIBRARY_PATH=/opt/fmrib/FreeSurfer_releases/6.0/lib/gcc/lib '...
        'mri_info ' anafile ' --vox2ras-tkr']);
    A=sscanf(result,'%f');
    Mvox2surf=reshape(A,4,4)'; % last couple of commands only got Mvox2Surf into Matlab format
    
    % Alon: Get transform conformed2World    
    [status,result] = system(['DYLD_LIBRARY_PATH=/opt/fmrib/FreeSurfer_releases/6.0/lib/gcc/lib '...
        'mri_info ' anafile ' --vox2ras']);
    A=sscanf(result,'%f');
    Mvox2space=reshape(A,4,4)';
    
    % get transform from RAS (surface) to World (mm, subject). 
    if (flipR & h==2)
        A=eye(4);
        A(1,1)=-1; A(1,4)=1;   % Flip is around +0.5 axis
        Msurf2space=Mvox2space*inv(Mvox2surf)*A;
    else
        Msurf2space=Mvox2space*inv(Mvox2surf);
    end;
    
    % Topology file
    % [vertex_coords, faces]
    [Vw,F]=read_surf([surf_dir filesep hem{h} '.white']);
    if (flipR & h==2)
        a=F(:,1);
        F(:,1)=F(:,3);
        F(:,3)=a;
    end;
    
    % Alon: I think the "closed" topology is only about the faces of the 
    % surface (i.e which vertices are connected). So this has nothing to do
    % with the coordinates of the vertices. 
    T1.perimeter_id={'CLOSED'};
    T1.encoding={'BINARY'};
    T1.data=F+1;                            % One-based indices
    toponame=[hem{h} '.CLOSED.topo'];
    caret_save([outd filesep toponame],T1);
    SPEC.data{end+1}=['CLOSEDtopo_file ' toponame];
       
    % Do all the Surfaces / coord files. 
    % Alon: this is about the vertex coords, not faces. Will save coord
    % files which store the coordinates of each vertex, but in subject
    % space (World, mm) coordinates rather than surface (RAS)
    for i=1:length(surfaceStr)
        filename=[surf_dir filesep hem{h} '.' surfaceStr{i}];
        if (exist(filename,'file'));
            [Vw]=read_surf(filename);
            C.configuration_id={surftype{i}};
            C.encoding={'BINARY'};
            C.data=[];
            if (noshift==1)
                C.data=Vw;
            else
                % C.data is nVertex x 3 (list of vertex coordinates in World (mm) coordinates)
                [C.data(:,1),C.data(:,2),C.data(:,3)]=spmj_affine_transform(Vw(:,1),Vw(:,2),Vw(:,3),Msurf2space);
            end;
            coord_name=[hem{h} '.' upper(surfaceStr{i}) '.coord'];
            caret_save([outd filesep coord_name],C);
            fprintf('wrote %s\n',coord_name);
            SPEC.data{end+1}=[surftype{i} 'coord_file ' coord_name];
        end;
    end;
    
    % Do the shape files and make then column is a surface_shape file
    if length(shape)>0
        S.num_cols=length(shape);
        S.encoding={'BINARY'};
        S.data=[];
        for i=1:length(shape)
            [S.data(:,i)]=shapesign(i).*read_curv([surf_dir filesep hem{h} '.' shape{i}]);
            S.column_name{i}=shape{i};
            S.column_color_mapping(i,:)=[-1.5000 1.5000];
        end;
        
        caret_save([outd filesep hem{h} '.surface_shape'],S);
        
        SPEC.data{end+1}=['surface_shape_file ' hem{h} '.surface_shape'];
    end;
    
    % Save the spec file
    specfilename=[outd filesep subj '.' hem{h} '.spec'];
    caret_save(specfilename,SPEC);
end;

setenv('SUBJECTS_DIR',old_dir);
