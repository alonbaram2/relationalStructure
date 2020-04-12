function T=caret_regionstats(regions,dependent,varargin);  
% T=caret_regionstats(regions,varargin);  
% INPUT:
%       regions: can be a) a paint file defining the regions 
%                b) a region structure with 
%                c) a matrix with one column per region of logical indices
%                d) a vector of node indices 
%       dependent: cell array of metric file names with 1 column per subj
% 
T=[];
P=caret_load(regions); 
reg=unique(P.data(:,1)); 
reg(reg==0)=[]; 

numdep=length(dependent);
for i=1:numdep; 
    Data{i}=caret_load(dependent{i}); 
    [dummy,name{i}]=fileparts(dependent{i});
    name{i}(name{i}=='.')='_';
end;

vararginoptions(varargin,{'name'});

numcols=size(Data{i}.data,2);

for i=1:numcols
    for r=1:length(reg) 
        D.SN=i; 
        D.region=reg(r); 
        for dep=1:numdep
            D.(name{dep})=nanmean(Data{dep}.data(P.data(:,1)==reg(r),i)); 
        end;
        T=addstruct(T,D);
    end;
end;    
    
varargout={T}; 