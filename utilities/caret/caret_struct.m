function M=caret_struct(kind,varargin);
% function M=caret_struct(kind,varargin);
%   Generates a new structure for caret-files
%   Examples:
%   M=caret_struct('coord','data',data);
%   T=caret_struct('metric','data',data);
%  INPUT: 
%       kind: 'coord','metric', RGB_paint','paint'
%  VARARIN: 
%     'data'
%     'num_nodes'
%     'num_cols'
%     'encoding'
%     'column_name'
%     'paintnames'
%     'scales'
%   ...
c=1;
M.encoding={'BINARY'};
while c<=length(varargin)
    switch varargin{c}
        case {'data','num_nodes','num_cols','encoding','column_name',...
                'paintnames','scales','configuration_id'}
            M=setfield(M,varargin{c},varargin{c+1});
            c=c+2;
        otherwise 
            error(sprintf('unknown tag %s',varargin{c})); 
    end;
end;

switch (kind)
    case 'coord'
        if (isfield(M,'data'))
            M.num_nodes=size(M.data,1);
            if size(M.data,2)~=3
                error('Coord-struct needs three columns');
            end;
        else
            error('Coord-struct needs data');
        end;
    case 'metric'
        if (isfield(M,'data'))
            M.num_rows=size(M.data,1);
            M.num_cols=size(M.data,2);
            if (~isfield(M,'column_name'))
                for c=1:M.num_cols
                    M.column_name{c}=sprintf('Column_%2.2d',c);
                end;
            end;
            if (~isfield(M,'column_color_mapping'))
                for c=1:M.num_cols
                    M.column_color_mapping(c,1:2)=[-5 5];
                end;
            end;
        else
            error('Metric-struct needs data');
        end;
    case 'RGB_paint'
        if (isfield(M,'data'))
            M.num_rows=size(M.data,1);
            M.num_cols=size(M.data,2)/3;
            if (~isfield(M,'column_name'))
                for c=1:M.num_cols
                    M.column_name{c}=sprintf('Column_%2.2d',c);
                end;
            end;
        end;
    case 'paint'
        if (~isfield(M,'data'))
            error('Paint-struct needs data');
        end;
        M.num_rows=size(M.data,1);
        M.num_cols=size(M.data,2);
        if (~isfield(M,'paintnames'))
            for c=0:(max(M.data(:)))
                M.paintnames{c+1}=sprintf('value_%2.2d',c);
            end;
        end;
        if (~isfield(M,'column_name'))
            for c=1:M.num_cols
                M.column_name{c}=sprintf('Column_%2.2d',c);
            end;
        end;
        M.num_paintnames=length(M.paintnames);
end;