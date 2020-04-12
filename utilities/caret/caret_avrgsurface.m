function caret_avrgsurface(varargin);
% function caret_avrgsurface(varargin);
% Makes an average surface structure for group statistics
% Assumes a set of Fiducial coord files that are aligned to a common
% topology, for example for the PALS atlas or the Free-surfer average
% Assumes that matlab is cd'd into caret_dir
%
% INPUT:
%   'avrgDir': Directory of the average surface: default 'fsaverage_sym'
%   'hemName': Names of the hemispheres (defaults: 'LeftHem','RightHem')
%   'hem':   : Naming convention for
%   'subjDir': Naming pattern of subject directory. Default ('xs*')
%   'hemis'  : hemispheres to do: default [1 2]
%   'coord'  : Name of the original coord files (e.g. 'WHITE', 'FIDUCIAL')  

% 24.09.2012 Tobi: added the fsaverage_sym option
% 20.01.2013 Joern: made universal - fsaverage or fsaverge_sym now signaled
%                   through the working directory and naming conventions 


avrgDir  = 'fsaverage_sym';
hemName  = {'LeftHem','RightHem'};
hem      = {'lh','rh'};
subjDir  = 'xs*';
topo_name ='';
coord     = 'WHITE';
hemis    = [1 2]; 

vararginoptions(varargin,{'avrgDir','topo_name','coord_name','subjDir','hemis'});

for h=1:2
    % Get average surfaces 
    dirMean=fullfile(avrgDir,hemName{h});
    topo_name=fullfile(dirMean,[hem{h} '.CLOSED.topo']);
    mean_coord=fullfile(dirMean,[hem{h} '.' coord '.coord']);
  
    % get individual surface banes 
    a=get_files(subjDir);
    for n=1:length(a)
        coord_name{n}=fullfile(a{n},hemName{h},[hem{h} '.' coord '.coord']);
    end;
    outname=fullfile(dirMean,[hem{h} '.avrgsurface.mat']);
    a=get_files('xs*');
    for n=1:length(a)
        coord_name{n}=fullfile(a{n},hemName{h},[hem{h} '.' coord '.coord']);
    end;
    outname=fullfile(dirMean,[hem{h} '.avrgsurface.mat']);
    
    % Calculate individual surfaces 
    N=length(coord_name);
    for i=1:N
        S{i}=caret_getsurface(coord_name{i},topo_name);
        S{i}=caret_calcarea(S{i});
        fprintf('Hem: %d Subject %d',h,i); 
    end;
    
    % Make the average surface
    MEAN=caret_getsurface(mean_coord,topo_name);
    MEAN=S{1};  
    
    % preallocate for speed
    P=size(S{1}.Nodes.area,1);          % Number of nodes 
    K=size(S{1}.Tiles.edgelength,1);    % Number of tiles 
    E=size(S{1}.Edges.Length,1);        % Number of Edges 
    NodeArea=zeros(P,N);
    TilesLength=zeros(K,3,N);
    EdgeLength=zeros(E,N);
    TileArea=zeros(K,N);

    % Calculates the indivudual areas and distances
    for i=1:N
        NodeArea(:,i)=S{i}.Nodes.area;
        TilesLength(:,:,i)=S{i}.Tiles.edgelength;
        EdgeLength(:,i)=S{i}.Edges.Length;
        TileArea(:,i)=S{i}.Tiles.area;
        fprintf('Hem: %d Subject %d',h,i); 
    end;
    
    % Average individual surfaces 
    MEAN.Nodes.area=mean(NodeArea,2);
    MEAN.Tiles.edgelength=mean(TilesLength,3);
    MEAN.Tiles.area=mean(TileArea,2);
    MEAN.Edges.Length=mean(EdgeLength,2);
    MEAN.A(3)=sum(MEAN.Nodes.area);
    MEAN.A(2)=sum(MEAN.Edges.Length(MEAN.Edges.num_occ==1))/2;
    save(outname,'-struct','MEAN');
end;