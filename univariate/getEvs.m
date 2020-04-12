function getEvs(rootData,subjects,iSub,glm)

% preprocessed data: should be 3D files (i.e. if using FSL for preprocessing as
% I did, you'll need to use fslsplit to split the 4D files to 3D files).
% This is so that SPM can deal with the files. An example 3D file name:
% /home/fs0/abaram/scripts/NN2020/preprocessed/sub-01/run1_cond1_0000.nii.
% As in the behavioural data, conditions (cond) 1&2 are -Corr blocks, while
% conditions 3&4 are +Corr blocks.
% My preprocessed images are unsmoothed, and I also registered each image
% to the downsampled (2x2x2mm) structural space of each subject. This step
% is not actually necessary - you can work in the functional space of each
% subject, though note some of the code will need adjusting.

mkdir(fullfile(rootData,'evs',glm,subjects{iSub}))
load([rootData '/bhv/processed/behMatlabVars.mat'],'rl','B','timings')
T = timings;

switch (glm)
    case('GLM1')
        nOnsetEvs = 6; % [AB],[C]@ both stim and outcome, [accept],[reject]@ choice
        for iRun=1:2
            for iCond=1:4
                names = cell(1,nOnsetEvs);
                onsets = cell(1,nOnsetEvs);
                durations = cell(1,nOnsetEvs);
                orth = cell(1,nOnsetEvs);
                
                % get variables of specific subject, will be easier to index
                % later
                
                subT.stim        = squeeze(T.stim(iSub,iRun,iCond,:));
                subT.outcome     = squeeze(T.outcome(iSub,iRun,iCond,:));
                subT.choice      = squeeze(T.choice(iSub,iRun,iCond,:));
                subRl.NAIVE.vCh  = squeeze(rl.NAIVE.vCh(iSub,iRun,iCond,:));
                subRl.STRUCT.vCh = squeeze(rl.STRUCT.vCh(iSub,iRun,iCond,:));
                subRt            = squeeze(B.rt(iSub,iRun,iCond,:));
                
                % regressors
                iEv = 1;
                
                names{iEv}     = 'AB_S';
                onsets{iEv}    = subT.stim(B.cues(iSub,iRun,iCond,:)~=2);
                
                % parametric modulators on stimulus onsets of AB
                pmod(iEv).name{1}  = 'vCh_NAIVE';
                pmod(iEv).param{1} = subRl.NAIVE.vCh(B.cues(iSub,iRun,iCond,:)~=2);
                pmod(iEv).param{1} = pmod(iEv).param{1} - mean(pmod(iEv).param{1}); % demean
                pmod(iEv).poly{1}  = 1;
                
                pmod(iEv).name{2}  = 'vCh_STRUCT';
                pmod(iEv).param{2} = subRl.STRUCT.vCh(B.cues(iSub,iRun,iCond,:)~=2);
                pmod(iEv).param{2} = pmod(iEv).param{2} - mean(pmod(iEv).param{2}); % demean
                pmod(iEv).poly{2}  = 1;
                
                % add reaction time as a neusance regressor
                pmod(iEv).name{3}  = 'rt';
                pmod(iEv).param{3} = subRt(B.cues(iSub,iRun,iCond,:)~=2);
                pmod(iEv).param{3} = pmod(iEv).param{3} - mean(pmod(iEv).param{3}); % demean
                pmod(iEv).poly{3}  = 1;
                
                iEv = iEv + 1;
                
                % other onset regressors
                names{iEv}     = 'C_S';
                onsets{iEv}    = subT.stim(B.cues(iSub,iRun,iCond,:)==2);
                iEv = iEv + 1;
                
                names{iEv}     = 'AB_O';
                onsets{iEv}    = subT.outcome(B.cues(iSub,iRun,iCond,:)~=2);
                iEv = iEv + 1;
                
                names{iEv}     = 'C_O';
                onsets{iEv}    = subT.outcome(B.cues(iSub,iRun,iCond,:)==2);
                iEv = iEv + 1;
                
                names{iEv}     = 'accept_Ch';
                onsets{iEv}    = subT.choice(B.choices(iSub,iRun,iCond,:)==1);
                iEv = iEv + 1;
                
                names{iEv}     = 'reject_Ch';
                onsets{iEv}    = subT.choice(B.choices(iSub,iRun,iCond,:)==0);
                iEv = iEv + 1;
                
                for iEv=1:length(names)
                    durations{iEv} = 0.1 .* ones(size(onsets{iEv}));  % all regressors are stick regressors (duration of 0.1 sec).
                    orth{iEv}      = false; % orthogonalise parametric modulators?
                end
                save(fullfile(rootData,'evs',glm,subjects{iSub},['run' int2str(iRun) '_cond' int2str(iCond)]),'names','onsets','durations','pmod','orth');
                
            end
        end
        
    case('GLM1_noRt')
        nOnsetEvs = 6; % [AB],[C]@ both stim and outcome, [accept],[reject]@ choice
        for iRun=1:2
            for iCond=1:4
                names = cell(1,nOnsetEvs);
                onsets = cell(1,nOnsetEvs);
                durations = cell(1,nOnsetEvs);
                orth = cell(1,nOnsetEvs);
                
                % get variables of specific subject, will be easier to index
                % later
                
                subT.stim        = squeeze(T.stim(iSub,iRun,iCond,:));
                subT.outcome     = squeeze(T.outcome(iSub,iRun,iCond,:));
                subT.choice      = squeeze(T.choice(iSub,iRun,iCond,:));
                subRl.NAIVE.vCh  = squeeze(rl.NAIVE.vCh(iSub,iRun,iCond,:));
                subRl.STRUCT.vCh = squeeze(rl.STRUCT.vCh(iSub,iRun,iCond,:));
                
                % regressors
                iEv = 1;
                
                names{iEv}     = 'AB_S';
                onsets{iEv}    = subT.stim(B.cues(iSub,iRun,iCond,:)~=2);
                
                % parametric modulators on stimulus onsets of AB
                pmod(iEv).name{1}  = 'vCh_NAIVE';
                pmod(iEv).param{1} = subRl.NAIVE.vCh(B.cues(iSub,iRun,iCond,:)~=2);
                pmod(iEv).param{1} = pmod(iEv).param{1} - mean(pmod(iEv).param{1}); % demean
                pmod(iEv).poly{1}  = 1;
                
                pmod(iEv).name{2}  = 'vCh_STRUCT';
                pmod(iEv).param{2} = subRl.STRUCT.vCh(B.cues(iSub,iRun,iCond,:)~=2);
                pmod(iEv).param{2} = pmod(iEv).param{2} - mean(pmod(iEv).param{2}); % demean
                pmod(iEv).poly{2}  = 1;
                
                iEv = iEv + 1;
                
                % other onset regressors
                names{iEv}     = 'C_S';
                onsets{iEv}    = subT.stim(B.cues(iSub,iRun,iCond,:)==2);
                iEv = iEv + 1;
                
                names{iEv}     = 'AB_O';
                onsets{iEv}    = subT.outcome(B.cues(iSub,iRun,iCond,:)~=2);
                iEv = iEv + 1;
                
                names{iEv}     = 'C_O';
                onsets{iEv}    = subT.outcome(B.cues(iSub,iRun,iCond,:)==2);
                iEv = iEv + 1;
                
                names{iEv}     = 'accept_Ch';
                onsets{iEv}    = subT.choice(B.choices(iSub,iRun,iCond,:)==1);
                iEv = iEv + 1;
                
                names{iEv}     = 'reject_Ch';
                onsets{iEv}    = subT.choice(B.choices(iSub,iRun,iCond,:)==0);
                iEv = iEv + 1;
                
                for iEv=1:length(names)
                    durations{iEv} = 0.1 .* ones(size(onsets{iEv}));  % all regressors are stick regressors (duration of 0.1 sec).
                    orth{iEv}      = false; % orthogonalise parametric modulators?
                end
                save(fullfile(rootData,'evs',glm,subjects{iSub},['run' int2str(iRun) '_cond' int2str(iCond)]),'names','onsets','durations','pmod','orth');
                
            end
        end
    case ('GLM2')
        nOnsetEvs = 8; % [A],[B],[C]@ both stim and outcome, [accept],[reject]@ choice
        for iRun=1:2
            for iCond=1:4
                names = cell(1,nOnsetEvs);
                onsets = cell(1,nOnsetEvs);
                durations = cell(1,nOnsetEvs);
                orth = cell(1,nOnsetEvs);
                
                % get variables of specific subject, will be easier to index
                % later
                
                subT.stim        = squeeze(T.stim(iSub,iRun,iCond,:));
                subT.outcome     = squeeze(T.outcome(iSub,iRun,iCond,:));
                subT.choice      = squeeze(T.choice(iSub,iRun,iCond,:));
                
                % regressors
                iEv = 1;
                
                names{iEv}     = 'A_S';
                onsets{iEv}    = subT.stim(B.cues(iSub,iRun,iCond,:)==0);
                iEv = iEv + 1;
                
                names{iEv}     = 'B_S';
                onsets{iEv}    = subT.stim(B.cues(iSub,iRun,iCond,:)==1);
                iEv = iEv + 1;
                
                names{iEv}     = 'C_S';
                onsets{iEv}    = subT.stim(B.cues(iSub,iRun,iCond,:)==2);
                iEv = iEv + 1;
                
                names{iEv}     = 'A_O';
                onsets{iEv}    = subT.outcome(B.cues(iSub,iRun,iCond,:)==0);
                iEv = iEv + 1;
                
                names{iEv}     = 'B_O';
                onsets{iEv}    = subT.outcome(B.cues(iSub,iRun,iCond,:)==1);
                iEv = iEv + 1;
                
                names{iEv}     = 'C_O';
                onsets{iEv}    = subT.outcome(B.cues(iSub,iRun,iCond,:)==2);
                iEv = iEv + 1;
                
                names{iEv}     = 'accept_Ch';
                onsets{iEv}    = subT.choice(B.choices(iSub,iRun,iCond,:)==1);
                iEv = iEv + 1;
                
                names{iEv}     = 'reject_Ch';
                onsets{iEv}    = subT.choice(B.choices(iSub,iRun,iCond,:)==0);
                iEv = iEv + 1;
                
                pmod = {};
                
                for iEv=1:length(names)
                    durations{iEv} = 0.1 .* ones(size(onsets{iEv}));  % all regressors are stick regressors (duration of 0.1 sec).
                    orth{iEv}      = false; % orthogonalise parametric modulators?
                end
                save(fullfile(rootData,'evs',glm,subjects{iSub},['run' int2str(iRun) '_cond' int2str(iCond)]),'names','onsets','durations','pmod','orth');
                
            end
        end
    case('GLM3')
        nOnsetEvs = 6; % [AB],[C]@ both stim and outcome, [accept],[reject]@ choice
        for iRun=1:2
            for iCond=1:4
                names = cell(1,nOnsetEvs);
                onsets = cell(1,nOnsetEvs);
                durations = cell(1,nOnsetEvs);
                orth = cell(1,nOnsetEvs);
                
                % get variables of specific subject, will be easier to index
                % later
                
                subT.stim         = squeeze(T.stim(iSub,iRun,iCond,:));
                subT.outcome      = squeeze(T.outcome(iSub,iRun,iCond,:));
                subT.choice       = squeeze(T.choice(iSub,iRun,iCond,:));
                subRl.STRUCT.crPe = squeeze(rl.STRUCT.crPe(iSub,iRun,iCond,:));
                
                % regressors
                iEv = 1;
                
                names{iEv}     = 'AB_O';
                onsets{iEv}    = subT.outcome(B.cues(iSub,iRun,iCond,:)~=2);
                
                % parametric modulator on outcome onsets of AB
                pmod(iEv).name{1}  = 'crPe_STRUCT';
                pmod(iEv).param{1} = subRl.STRUCT.crPe(B.cues(iSub,iRun,iCond,:)~=2);
                pmod(iEv).param{1} = pmod(iEv).param{1} - mean(pmod(iEv).param{1}); % demean
                pmod(iEv).poly{1}  = 1;
                
                iEv = iEv + 1;
                
                names{iEv}     = 'C_O';
                onsets{iEv}    = subT.outcome(B.cues(iSub,iRun,iCond,:)==2);
                iEv = iEv + 1;
                
                names{iEv}     = 'AB_S';
                onsets{iEv}    = subT.stim(B.cues(iSub,iRun,iCond,:)~=2);
                iEv = iEv + 1;
                
                names{iEv}     = 'C_S';
                onsets{iEv}    = subT.stim(B.cues(iSub,iRun,iCond,:)==2);
                iEv = iEv + 1;
                
                names{iEv}     = 'accept_Ch';
                onsets{iEv}    = subT.choice(B.choices(iSub,iRun,iCond,:)==1);
                iEv = iEv + 1;
                
                names{iEv}     = 'reject_Ch';
                onsets{iEv}    = subT.choice(B.choices(iSub,iRun,iCond,:)==0);
                iEv = iEv + 1;
                
                for iEv=1:length(names)
                    durations{iEv} = 0.1 .* ones(size(onsets{iEv}));  % all regressors are stick regressors (duration of 0.1 sec).
                    orth{iEv}      = false; % orthogonalise parametric modulators?
                end
                save(fullfile(rootData,'evs',glm,subjects{iSub},['run' int2str(iRun) '_cond' int2str(iCond)]),'names','onsets','durations','pmod','orth');
                
            end
        end
        
end
