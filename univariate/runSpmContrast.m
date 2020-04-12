function runSpmContrast(rootData,sub,glm)
% # calculate contrasts
cd(fullfile(rootData, 'spmDirs',glm,sub));
load SPM
% delete preivious calculations of contrasts, if they exist
if isfield(SPM,'xCon')
    SPM=rmfield(SPM,'xCon');
end
delete('*spmT*')
delete('con*')


switch(glm)
    case 'GLM1'
        % contrasts names
        cNames = {'STRUCT','STRUCT>NAIVE','STRUCT_neg','NAIVE','rt'};
        con=[];
        
        % contrast 1: STRUCT
        c=zeros(1,size(SPM.xX.X,2));
        for r=1:8 % blocks
            v = SPM.Sess(r).col; % all indeces of regressors of block
            c(v(3))=1; % AB_S x vCh_STRUCT
        end
        con = [con; c];
        
        % contrast 2: STRUCT>NAIVE
        c=zeros(1,size(SPM.xX.X,2));
        for r=1:8
            v = SPM.Sess(r).col;
            c(v(3))=1; % AB_S x vCh_STRUCT
            c(v(2))=-1; % AB_S x vCh_NAIVE
        end
        con = [con; c];
        
        % contrast 3: -STRUCT
        c=zeros(1,size(SPM.xX.X,2));
        for r=1:8
            v = SPM.Sess(r).col;
            c(v(3))=-1; % - AB_S x vCh_STRUCT
        end
        con = [con; c];
        
        % contrast 4: NAIVE
        c=zeros(1,size(SPM.xX.X,2));
        for r=1:8
            v = SPM.Sess(r).col;
            c(v(2))= 1; % AB_S x vCh_NAIVE
        end
        con = [con; c];
        
        % contrast 5: rt
        c=zeros(1,size(SPM.xX.X,2));
        for r=1:8
            v = SPM.Sess(r).col;
            c(v(4))=1; % AB_S x rt
        end
        con = [con; c];
        
    case 'GLM1_noRt'
        % contrasts names
        cNames = {'STRUCT','STRUCT>NAIVE','STRUCT_neg','NAIVE'};
        con=[];
        
        % contrast 1: STRUCT
        c=zeros(1,size(SPM.xX.X,2));
        for r=1:8 % blocks
            v = SPM.Sess(r).col; % all indeces of regressors of block
            c(v(3))=1; % AB_S x vCh_STRUCT
        end
        con = [con; c];
        
        % contrast 2: STRUCT>NAIVE
        c=zeros(1,size(SPM.xX.X,2));
        for r=1:8
            v = SPM.Sess(r).col;
            c(v(3))=1; % AB_S x vCh_STRUCT
            c(v(2))=-1; % AB_S x vCh_NAIVE
        end
        con = [con; c];
        
        % contrast 3: -STRUCT
        c=zeros(1,size(SPM.xX.X,2));
        for r=1:8
            v = SPM.Sess(r).col;
            c(v(3))=-1; % - AB_S x vCh_STRUCT
        end
        con = [con; c];
        
        % contrast 4: NAIVE
        c=zeros(1,size(SPM.xX.X,2));
        for r=1:8
            v = SPM.Sess(r).col;
            c(v(2))= 1; % AB_S x vCh_NAIVE
        end
        con = [con; c];
        
    case 'GLM3'
        % contrasts names
        cNames = {'crPe_STRUCT'};
        con=[];
        
        % contrast 1: STRUCT
        c=zeros(1,size(SPM.xX.X,2));
        for r=1:8 % blocks
            v = SPM.Sess(r).col; % all indeces of regressors of block
            c(v(2))=1; % AB_S x vCh_STRUCT
        end
        con = [con; c];        
               
end

% caclulate contrast images
for iCon=1:numel(cNames)
    SPM.xCon(iCon)=spm_FcUtil('Set',sprintf('%s',cNames{iCon}), 'T', 'c',con(iCon,:)',SPM.xX.xKXs);
end
SPM=spm_contrasts(SPM,[1:length(SPM.xCon)]);
save SPM SPM;
