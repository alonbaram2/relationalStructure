function S=caret_makereferencestrcut(varargin)
% function S=spmj_makereferencestrcut(varargin)
% based on a location structure.
P=1;
s=1;
while (~isempty(P))
    P=spm_get(inf,'*.metric',sprintf('Get metric files for subject %d',s));
    if (~isempty(P))
        S.images{s}=P;
        s=s+1;
    end;
end;
S.deformations=spm_get(inf,'d*.coord',sprintf('Get deformed coordinate files'));
varargout={S};
