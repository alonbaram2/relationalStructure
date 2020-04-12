function T=caret_regionstat(varargin);
%  function T=caret_regionstat(varargin);
%       Extracts the values of images/metric over the regions or points 
% 'images': A cell array of images / metric files 
% 'deformations': Deformations from atlas space into original space. For
%           caret, this can be deformed_fiducial metric file, if not given,
%           no deformation is assumed 
% 'regions': cell array of region-structures 
%           R.type='node';
%           R.location=[x,y,z];
%           --------------------
%           R.type='node'
%           R.location=[x,y,z];
%           -------------------- 
%           R.type='paint'
%           R.location=[filename];
% -------------------------------------------------------------------------
% v.1.0 Joern Diedrichsen jdiedric@bme.jhu.edu
% v.1.1 Made much more universal: points are possible, reference back to
%   the original space are possible and data extraction of the full time
%   series 
% 

% Result structure 
T.region=[];
T.region_name={};
T.image=[];
T.image_name={};
T.column=[];
T.column_name=[];
T.reg_size=[];
T.data=[];

images=[];
regions=[];
deformations=[];

vararginoptions(varargin,{'images','regions','deformations'});

% -------------------------------------
% Prompt with GUI if arguments are empty 
if (isempty(regions)) 
    regions=caret_getregions();
end;
if (isempty(images))
    images=spm_get(inf,'.metric',sprintf('Select images / metric files'));
end;
if (iscell(images))
    images=char(images);
end;
% ---------------------------------------
% Calculate the regions 
LOC=caret_calcregions(regions);

% ---------------------------------------
% Loop images and columns 
num_regions=length(LOC);
for i=1:size(images,1)
    M=caret_load(deblank(images(i,:)));
    fprintf('Metric: %s\n',deblank(images(i,:)));
    for r=1:num_regions
        for c=1:M.num_cols
            T.region(end+1,1)=r;
            T.image(end+1,1)=i;
            T.image_name{end+1,1}=images(i,:);
            T.column(end+1,1)=c;
            T.column_name{end+1,1}=M.column_name{c};
            T.reg_size(end+1,1)=length(LOC{r});
            d=M.data(LOC{r},c); 
            if sum(~isnan(d))==0
                T.data(end+1,1)=NaN;
            else 
                T.data(end+1,1)=nanmean(d);
            end;
        end;    
    end;
end;
