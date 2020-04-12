function LOC=spmj_calcregions(R,varargin)
% function LOC=spmj_calcregions(R,varargin)
% returns a cell array of coordinate locations to be sampled
% based on a location structure.
gridwidth=2;        % 2mm grid width for sample_image
if (~iscell(R))
    S{1}=R;
    R=S;
end;
LOC={};
for c=1:length(R)
    switch  (R{c}.type)
        case 'point'
            for i=1:size(R{c}.location,1);
                LOC{end+1}=R{c}.location(i,:);
            end;
        case 'sphere'
            % make sphere
            left=-R{c}.radius.*[1 1 1];
            right=+R{c}.radius.*[1 1 1];
            s=R{c}.radius/(ceil(R{c}.radius/gridwidth));
            [x,y,z]=meshgrid(left(1):s:right(1),left(2):s:right(2),left(3):s:right(3));
            x=x(:);y=y(:);z=z(:);
            rad=sqrt(x.^2+y.^2+z.^2);
            indx=find(rad<=R{c}.radius);
            x=x(indx);y=y(indx);z=z(indx);
            % add for all locations
            for i=1:size(R{c}.location,1);
                cent=R{c}.location(i,:);
                LOC{end+1}=[x+cent(1) y+cent(2) z+cent(3)];
            end;
        case 'box'
            % make Box
            left=-R{c}.size/2;
            right=+R{c}.size/2;
            [x,y,z]=meshgrid(left(1):gridwidth:right(1),left(2):gridwidth:right(2),left(3):gridwidth:right(3));
            x=x(:);y=y(:);z=z(:);
            for i=1:size(R{c}.location,1);
                cent=R{c}.location(i,:);
                LOC{end+1}=[x+cent(1) y+cent(2) z+cent(3)];
            end;
        case 'image'
            % Load images
            V=spm_vol(R{c}.location);
            for v=1:length(V)
                X=spm_read_vols(V(v));
                [x,y,z]=ind2sub(size(X),find(X>0));
                [x,y,z]=spmj_affine_transform(x,y,z,V(v).mat);
                LOC{end+1}=[x+cent(1) y+cent(2) z+cent(3)];
            end;
        case 'paint'
        case 'node'
            LOC{end+1}=R{c}.location;
    end;
end;