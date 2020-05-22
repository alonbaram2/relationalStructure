clear
rootScripts = '/home/fs0/abaram/scripts/NN2020';
rootData    = '/home/fs0/abaram/scratch/NN2020';
subjects = {'sub-01','sub-02','sub-03','sub-04','sub-05','sub-06','sub-07','sub-08','sub-09',...
    'sub-10','sub-11','sub-12','sub-13','sub-14','sub-15','sub-16','sub-17','sub-18','sub-19',...
    'sub-20','sub-21','sub-22','sub-23','sub-24','sub-25','sub-26','sub-27','sub-28'};

nSub = length(subjects);
nRun  = 2; % two independent runs, each with the 4 conditions
nCond = 4; % the 4 experimental conditions
nTrialsPerBlock = 30;
nTrialsPerBlock_instruct = 42;

%% create behaviour MATLAB structure (B)
% B dims are 28x2x4x30: subject x runs x conditions x number of trials per scanning block

B.cues     = nan(nSub,nRun,nCond,nTrialsPerBlock); 
B.outcomes = nan(nSub,nRun,nCond,nTrialsPerBlock); 
B.choices  = nan(nSub,nRun,nCond,nTrialsPerBlock);
B.probs    = nan(nSub,nRun,nCond,nTrialsPerBlock); 
B.rt       = nan(nSub,nRun,nCond,nTrialsPerBlock); 

timings.stim    = nan(nSub,nRun,nCond,nTrialsPerBlock);
timings.outcome = nan(nSub,nRun,nCond,nTrialsPerBlock);
timings.choice  = nan(nSub,nRun,nCond,nTrialsPerBlock);

B_ins.cues     = nan(nSub,nRun,nCond,nTrialsPerBlock_instruct); 
B_ins.outcomes = nan(nSub,nRun,nCond,nTrialsPerBlock_instruct); 
B_ins.choices  = nan(nSub,nRun,nCond,nTrialsPerBlock_instruct);
B_ins.rt = nan(nSub,nRun,nCond,nTrialsPerBlock_instruct); 

for iSub = 1:length(subjects)
    sub = subjects{iSub};
    f = [rootData '/bhv/raw/scanDay/' sub '/scan.mat'];
    load(f);
    for run = 1:nRun % two runs, each with the 4 conditions
        for cond = 1:nCond % the 4 experimental conditions, already reordered:
                       % 1: stimSet 1, -Corr.
                       % 2: stimSet 2, -Corr. 
                       % 3: stimSet 1, +Corr. 
                       % 4: stimSet 2, +Corr.
                       
            B.cues(iSub,run,cond,:)     = data.run{run}.cond{cond}.cues;    % 0: A. 1: B. 2: C (A & B are related. C is control). 
            B.outcomes(iSub,run,cond,:) = data.run{run}.cond{cond}.rewards; % 1: king of hearts. 0: two of spades. 0 will change to -1 later to comply with the equations of the RL model. 
            B.choices(iSub,run,cond,:)  = data.run{run}.cond{cond}.choices; % 1: accept. 0: reject. 
            B.probs(iSub,run,cond,:)    = data.run{run}.cond{cond}.probs; % congruency between choice and outcome.            
            B.rt(iSub,run,cond,:)       = data.run{run}.cond{cond}.RT;      % reaction times. 
                        
            timings.stim(iSub,run,cond,:)    = data.run{run}.cond{cond}.stimOn-data.run{run}.cond{cond}.tTrig1;    % stimulus presentation onset
            timings.outcome(iSub,run,cond,:) = data.run{run}.cond{cond}.feedbackOn-data.run{run}.cond{cond}.tTrig1; % outcome presentation onset
            timings.choice(iSub,run,cond,:)  = data.run{run}.cond{cond}.tKeyPress-data.run{run}.cond{cond}.tTrig1; % choice/button press time
        end
    end
    
    f = [rootData '/bhv/raw/scanDay/' sub '/instruct.mat'];
    load(f);
    run = 1; % only one run in instruction session
    for cond = 1:nCond % the 4 experimental conditions, already reordered:
        % 1: stimSet 1, -Corr.
        % 2: stimSet 2, -Corr.
        % 3: stimSet 1, +Corr.
        % 4: stimSet 2, +Corr.
        
        B_ins.cues(iSub,run,cond,:)     = data.run{run}.cond{cond}.cues;    % 0: A. 1: B. 2: C (A & B are related. C is control).
        B_ins.outcomes(iSub,run,cond,:) = data.run{run}.cond{cond}.rewards; % 1: king of hearts. 0: two of spades
        B_ins.choices(iSub,run,cond,:)  = data.run{run}.cond{cond}.choices; % 1: accept. 0: reject.
        B_ins.rt(iSub,run,cond,:)       = data.run{run}.cond{cond}.RT;      % reaction times.
    end    
end

% change negative outcomes from 0 to -1 to fit with the assumptions of the
% RL model.
B.outcomes(B.outcomes==0)         = -1; 
B_ins.outcomes(B_ins.outcomes==0) = -1; 

save([rootData '/bhv/processed/behMatlabVars.mat'],'B','B_ins','timings');

%% Fit and run RL models
% To increase power in the fitting process, we first concatenate all data
% from blocks with the same correlation (-/+Corr) from the prescanning and
% scanning sessions.
% These are not the cross-validated models that produce figures 2a and 2b
% (these wil;l comew later). THe models fitted here will be the ones used
% forthe fMRI analyses. 

R = 20; % number of random initialisations of parameters, to ensure fitting is not local minima
theta_NAIVE  = nan(nSub,2,2); % 2 block types (+/-Corr, collapse over different stimuli sets). 2 parameters to fit (alpha, beta) 
theta_STRUCT = nan(nSub,2,5); % 2 block types (+/-Corr, collapse over different stimuli sets). 5 parameters to fit (alpha, beta, 3 cross-terms) 
negLogLiklihood_NAIVE = nan(nSub,2); % 2 block types (+/-Corr, collapse over different stimuli sets). 
negLogLiklihood_STRUCT = nan(nSub,2); % 2 block types (+/-Corr, collapse over different stimuli sets).

for iSub = 1:length(subjects)
    sub = subjects{iSub};

    % Get a vector of trial numbers of last trial of blocks - model
    % estimations should reset at the beginning of blocks. should be
    % [42, 84, 114, 144, 174, 204]
    blocks = [nTrialsPerBlock_instruct,2*nTrialsPerBlock_instruct,...
        2*nTrialsPerBlock_instruct + nTrialsPerBlock:nTrialsPerBlock:2*nTrialsPerBlock_instruct + 4*nTrialsPerBlock];
    
    for iCond=[1,3]; % loop over neg corr and then pos corr    
        % concatenate data from blocks with same correlation (iCond and
        % iCond+1: blocks 1,2 are -Corr, blocks 3,4 are +Corr)
        cues_cat     = [squeeze(B_ins.cues(iSub,1,iCond,:))',squeeze(B_ins.cues(iSub,1,iCond+1,:))',squeeze(B.cues(iSub,1,iCond,:))',squeeze(B.cues(iSub,1,iCond+1,:))',squeeze(B.cues(iSub,2,iCond,:))',squeeze(B.cues(iSub,2,iCond+1,:))'];
        outcomes_cat = [squeeze(B_ins.outcomes(iSub,1,iCond,:))',squeeze(B_ins.outcomes(iSub,1,iCond+1,:))',squeeze(B.outcomes(iSub,1,iCond,:))',squeeze(B.outcomes(iSub,1,iCond+1,:))',squeeze(B.outcomes(iSub,2,iCond,:))',squeeze(B.outcomes(iSub,2,iCond+1,:))'];
        choices_cat  = [squeeze(B_ins.choices(iSub,1,iCond,:))',squeeze(B_ins.choices(iSub,1,iCond+1,:))',squeeze(B.choices(iSub,1,iCond,:))',squeeze(B.choices(iSub,1,iCond+1,:))',squeeze(B.choices(iSub,2,iCond,:))',squeeze(B.choices(iSub,2,iCond+1,:))'];
        
        % store maximum likelihood parameters for NAIVE: [alpha,beta].
        [theta_NAIVE(iSub,(iCond+1)/2,:), negLogLiklihood_NAIVE(iSub,(iCond+1)/2)] = fitRlParams(cues_cat,outcomes_cat,choices_cat,false,blocks,R);
        % store maximum likelihood parameters for STRUCT: [alpha,beta,H_AB,H_AC,H_BC].
        [theta_STRUCT(iSub,(iCond+1)/2,:), negLogLiklihood_STRUCT(iSub,(iCond+1)/2)] = fitRlParams(cues_cat,outcomes_cat,choices_cat,true,blocks,R);
    end
end

thetaMl.NAIVE    = theta_NAIVE;
thetaMl.STRUCT   = theta_STRUCT;
negLogLik.NAIVE  = negLogLiklihood_NAIVE;
negLogLik.STRUCT = negLogLiklihood_STRUCT;
negLogLik.N      = numel(cues_cat)* 2 * nSub; % total number of examples (trials), will be used to calculate AIC and BIC

% save variables
save([rootData '/bhv/processed/behMatlabVars.mat'],'thetaMl','negLogLik','-append');



%% get trial-by-trial values of RL variables using STRUCT and NAIVE models 
%with subject-specific fitted parameters

% RL variables that are tracked (for both NAIVE and STRUCT models): 

% g      : Outcome probability estimation.
% vCh    : Value of the chosen action, which always has the same magnitude 
%          but has an opposite sign to g on rejected trials
% epsilon: Outcome prediction error
% crPe   : The "correctness" prediction error, which has the same magnitude
%           as epsilon but a sign that depends on the congruence of the 
%           outcome and the subject's choice.
load([rootData '/bhv/processed/behMatlabVars.mat'],'thetaMl','negLogLik','B','B_ins')

rl.NAIVE.g       = nan(nSub,nRun,nCond,nTrialsPerBlock);
rl.NAIVE.vCh     = nan(nSub,nRun,nCond,nTrialsPerBlock);
rl.NAIVE.epsilon = nan(nSub,nRun,nCond,nTrialsPerBlock);
rl.NAIVE.crPe    = nan(nSub,nRun,nCond,nTrialsPerBlock);

rl.STRUCT.g       = nan(nSub,nRun,nCond,nTrialsPerBlock);
rl.STRUCT.vCh     = nan(nSub,nRun,nCond,nTrialsPerBlock);
rl.STRUCT.epsilon = nan(nSub,nRun,nCond,nTrialsPerBlock);
rl.STRUCT.crPe    = nan(nSub,nRun,nCond,nTrialsPerBlock);

% cues are coded as 0,1,2; for Matlab indexing move to be 1,2,3
cues = B.cues+1;

for iSub=1:nSub
    for iRun=1:nRun
        for iCond=1:nCond
            
            alpha_NAIVE  = squeeze(thetaMl.NAIVE(iSub,ceil(iCond/2),1));
            alpha_STRUCT = squeeze(thetaMl.STRUCT(iSub,ceil(iCond/2),1));
            crossTerms   = squeeze(thetaMl.STRUCT(iSub,ceil(iCond/2),3:5))';
            
            % initialise g at 0 for all three stimuli
            g_NAIVE = zeros(3,1); % 3 stimuli
            g_STRUCT = zeros(3,1); % 3 stimuli
            
            for t=1:nTrialsPerBlock
                % find current cue. Cues are coded as 0,1,2, so need to add
                % +1 for Matlab indexing               
                c = B.cues(iSub,iRun,iCond,t) +1;
                
                % calculate outcome probability etimation and outcome
                % prediciton error
                rl.NAIVE.g(iSub,iRun,iCond,t)        = g_NAIVE(c);
                rl.NAIVE.epsilon(iSub,iRun,iCond,t)  = B.outcomes(iSub,iRun,iCond,t) - g_NAIVE(c);
                rl.STRUCT.g(iSub,iRun,iCond,t)       = g_STRUCT(c);
                rl.STRUCT.epsilon(iSub,iRun,iCond,t) = B.outcomes(iSub,iRun,iCond,t) - g_STRUCT(c);                
                
                % flip sign if subject rejected for vCh and crPe.
                if B.choices(iSub,iRun,iCond,t)==1 % trial accepted
                    rl.NAIVE.vCh(iSub,iRun,iCond,t)   = rl.NAIVE.g(iSub,iRun,iCond,t);
                    rl.NAIVE.crPe(iSub,iRun,iCond,t)  = rl.NAIVE.epsilon(iSub,iRun,iCond,t);
                    rl.STRUCT.vCh(iSub,iRun,iCond,t)  = rl.STRUCT.g(iSub,iRun,iCond,t);
                    rl.STRUCT.crPe(iSub,iRun,iCond,t) = rl.STRUCT.epsilon(iSub,iRun,iCond,t);                      
                elseif B.choices(iSub,iRun,iCond,t)==0 % trial rejected
                    rl.NAIVE.vCh(iSub,iRun,iCond,t)   = - rl.NAIVE.g(iSub,iRun,iCond,t);
                    rl.NAIVE.crPe(iSub,iRun,iCond,t)  = - rl.NAIVE.epsilon(iSub,iRun,iCond,t);                    
                    rl.STRUCT.vCh(iSub,iRun,iCond,t)  = - rl.STRUCT.g(iSub,iRun,iCond,t);
                    rl.STRUCT.crPe(iSub,iRun,iCond,t) = - rl.STRUCT.epsilon(iSub,iRun,iCond,t);                    
                end
                g_NAIVE  = getRlValue(g_NAIVE,B.outcomes(iSub,iRun,iCond,t),B.cues(iSub,iRun,iCond,t),3,alpha_NAIVE,[0,0,0]);
                g_STRUCT = getRlValue(g_STRUCT,B.outcomes(iSub,iRun,iCond,t),B.cues(iSub,iRun,iCond,t),3,alpha_STRUCT,crossTerms);

            end
        end
    end
end

% save
save([rootData '/bhv/processed/behMatlabVars.mat'],'rl','-append');


%% stats and plots
%
% Subjects use the STRUCT model: plots and formal model comparison
%

load([rootData '/bhv/processed/behMatlabVars.mat'],'thetaMl','negLogLik','B','rl','B_ins')

%###
% fitted cross terms plot - fig *****

figure
subplot(2,1,1) % plot +Corr blocks cross-terms
hold on
title('+Corr blocks')
h1 = boxplot(squeeze(thetaMl.STRUCT(:,2,3:5)));
set(h1,'LineWidth',2)
set(gca,'YTick',[0,1],'FontSize',16,'XTick',[])
plot([0.5,3.5],[0,0],'--k','Linewidth',1.5)
plot([0.5,3.5],[1,1],'-k','Linewidth',1)
ylim([-0.5,1.3])

subplot(2,1,2) % plot -Corr blocks cross-terms
hold on
title('-Corr blocks`')
h2 = boxplot(squeeze(thetaMl.STRUCT(:,1,3:5)));
set(h2,'LineWidth',2)
set(gca,'YTick',[-1,0],'FontSize',16,'XTickLabels',{'H_{AB}','H_{AC}','H_{BC}'},'TickLabelInterpreter','tex')
plot([0.5,3.5],[0,0],'--k','Linewidth',1.5)
plot([0.5,3.5],[-1,-1],'-k','Linewidth',1)
ylim([-1.3,0.5])

% fitted cross-terms (STRUCT) stats
[~,P,~,S] = ttest(abs(thetaMl.STRUCT(:,1,3)) + abs(thetaMl.STRUCT(:,2,3)) , ...
    (abs(thetaMl.STRUCT(:,1,4))+abs(thetaMl.STRUCT(:,1,5)))./2 + (abs(thetaMl.STRUCT(:,2,4))+abs(thetaMl.STRUCT(:,2,5)))./2 , 'tail','right');
stats.crossTerms.P = P;
stats.crossTerms.t = S.tstat;
fprintf('Test cross-terms Corr vs 0Corr: P=%g, t=%g\n',P,S.tstat);

% ###
% model comparison
N = negLogLik.N; % number of examples (trials)

disp('model comparison: smaller is better')
% NAIVE
LL = sum(sum( - negLogLik.NAIVE)); % LL: log likelihood. collapse over -/+ Corr blocks and over subjects
k = 2 * 2 * nSub; %number of parameters fitted (2 per subject: alpha and beta) * 2 structures * nSub
AIC = 2*k - 2*LL;
BIC = log(N)*k - 2*LL;
disp(sprintf('NAIVE: -LL: %g , AIC: %g ,  BIC: %g\n' ,-LL,AIC,BIC))
stats.modelComparison.NAIVE.negLL = -LL;
stats.modelComparison.NAIVE.AIC = AIC;
stats.modelComparison.NAIVE.BIC = BIC;

% STRUCT
LL = sum(sum( - negLogLik.STRUCT)); % LL: log likelihood. collapse over -/+ Corr blocks and over subjects
k = 5 * 2 * nSub; %number of parameters fitted (5 per subject: alpha, beta and 3 cross terms) * 2 structures * nSub
AIC = 2*k - 2*LL;
BIC = log(N)*k - 2*LL;
disp(sprintf('STRUCT: -LL: %g , AIC: %g ,  BIC: %g\n' ,-LL,AIC,BIC))
stats.modelComparison.STRUCT.negLL = -LL;
stats.modelComparison.STRUCT.AIC = AIC;
stats.modelComparison.STRUCT.BIC = BIC;

% #############
%
%% Possible confounds: accuracy, reaction time and switch costs analyses
%

load([rootData '/bhv/processed/behMatlabVars.mat'],'B','rl')

% 5 regs of interest (relationalStructure, pe x structure interaction, latentVariable (AB vs C)
% vCh_STRUCT(value of chosen action), vCh_DUFF (vCh_STRUCT - vCh_Naive_,
% 6 confound regressors (rt, log(rt), correct, switch cost x 3)
stats.confoundsCorr.corrCoefs = nan(nSub,5,6);
stats.confoundsCorr.interestRegsNames = {'relationalStruct','peXstruct','latentVariable (AB v C)','vCh STRUCT','vCh DIFF'};
stats.confoundsCorr.confoundRegsNames = {'rt', 'log(rt)', 'correct', 'taskSwitch1','taskSwitch2','taskSwitch3'};

stats.accuracy.means = nan (nSub,3); % percent correct for 3 conditions: [+,-,0] Corr
stats.logrt.means = nan (nSub,3); % mean log(rt) for 3 conditions: [+,-,0] Corr
stats.logrtByTaskSwitchConds.means = nan (nSub,3); % mean log(rt) for the task switching conditions: same: A=>A, B=>B, C=>C; related: A=>B, B=>A; control: A=>C, B=>C, C=>A, C=>B

% save confound values to use later to construct confound EVs for GLMs. 
confounds.rt          = nan(size(B.cues)); 
confounds.logrt       = nan(size(B.cues)); 
confounds.correct     = nan(size(B.cues)); 
confounds.switchTask1 = nan(size(B.cues)); 
confounds.switchTask2 = nan(size(B.cues)); 
confounds.switchTask3 = nan(size(B.cues)); 

for iSub=1:length(subjects)
    
    % find indeces of different trial types: 
    
    % First find indeces of all -Corr blocks:    
    negCorrBlocksInds = false(size(B.cues(iSub,:,:,:)));
    negCorrBlocksInds(1,:,1:2,:) = true;
    % find all A and B cues in -Corr blocks
    negCorrTrialsInds = negCorrBlocksInds & (B.cues(iSub,:,:,:) ~= 2 );            
    % find indeces of all +Corr blocks:    
    posCorrBlocksInds = false(size(B.cues(iSub,:,:,:)));
    posCorrBlocksInds(1,:,3:4,:) = true;
    % find all A and B cues in +Corr blocks
    posCorrTrialsInds = posCorrBlocksInds & (B.cues(iSub,:,:,:) ~= 2 );          
    % find all C cues
    zeroCorrTrialsInds = false(size(B.cues(iSub,:,:,:)));
    zeroCorrTrialsInds(B.cues(iSub,:,:,:)==2) = true;
        
    % find the correct choice trials according to STRUCT model    
    correct =  B.choices(iSub,:,:,:) == round(B.probs(iSub,:,:,:)); % change choices to be {-1,1} instead of {0,1}
    
    subjectPeSTRUCT = rl.STRUCT.crPe(iSub,:,:,:);
    
    % store rt of specific subject in a temporary variable 
    subjectRts = B.rt(iSub,:,:,:);      
    % percent correct in +Corr trials
    stats.accuracy.means(iSub,1) = sum(correct(posCorrTrialsInds)) ./ length(correct(posCorrTrialsInds));
    % mean log(rt) in +Corr
    stats.logrt.means(iSub,1) = mean(log(subjectRts(posCorrTrialsInds)));       
    % percent correct in -Corr trials
    stats.accuracy.means(iSub,2) = sum(correct(negCorrTrialsInds)) ./ length(correct(negCorrTrialsInds));
    % mean log(rt) in -Corr
    stats.logrt.means(iSub,2) = mean(log(subjectRts(negCorrTrialsInds)));    
    % percent correct in 0Corr trials
    stats.accuracy.means(iSub,3) = sum(correct(zeroCorrTrialsInds)) ./ length(correct(zeroCorrTrialsInds));
    % mean log(rt) in 0Corr
    stats.logrt.means(iSub,3) = mean(log(subjectRts(zeroCorrTrialsInds)));
    
    % split trials depending on the preceding trial: 0: same (A=>A, B=>B,
    % C=>C); 1: related (A=>B, B=>A); control (A=>C,B=>C,C=>A,C=>B).     
    switchTask1 = nan(size(subjectRts));
    for iRun=1:size(switchTask1,2)
        for iCond = 1:size(switchTask1,3)
            for iTrial = 2:size(switchTask1,4)                
                if B.cues(iSub,iRun,iCond,iTrial) == B.cues(iSub,iRun,iCond,iTrial-1)
                    switchTask1(1,iRun,iCond,iTrial) = 0; % same  
                elseif  (B.cues(iSub,iRun,iCond,iTrial) == 0 & B.cues(iSub,iRun,iCond,iTrial-1) ==1) | ...
                        (B.cues(iSub,iRun,iCond,iTrial) == 1 & B.cues(iSub,iRun,iCond,iTrial-1) ==0)
                    switchTask1(1,iRun,iCond,iTrial) = 1; % related 
                else
                    switchTask1(1,iRun,iCond,iTrial) = 2; % control
                end                    
            end
        end
    end
    % Store mean of of log(rt), arranged by the task switching condition
    % (same, related, control). 
    stats.logrtByTaskSwitchConds.means(iSub,1)    = mean(log(subjectRts(switchTask1==0)));
    stats.logrtByTaskSwitchConds.means(iSub,2) = mean(log(subjectRts(switchTask1==1)));
    stats.logrtByTaskSwitchConds.means(iSub,3) = mean(log(subjectRts(switchTask1==2)));    
    
    % get task switch confound regressors: check all 3 valid options 
    % (same<related<control - thats switchTask1; [same,related]<control; same<[related,control]
    switchTask2=nan(size(switchTask1));
    switchTask2(switchTask1==0 | switchTask1==1) = -1;
    switchTask2(switchTask1==2) = 1; 
    switchTask3=nan(size(switchTask1));
    switchTask3(switchTask1==0) = -1;
    switchTask3(switchTask1==1 | switchTask1==2) = 1; 
         
    confounds.rt(iSub,:,:,:)          = subjectRts;
    confounds.logrt(iSub,:,:,:)       = log(subjectRts);
    confounds.correct(iSub,:,:,:)     = correct;
    confounds.switchTask1(iSub,:,:,:) = switchTask1;
    confounds.switchTask2(iSub,:,:,:) = switchTask2;
    confounds.switchTask3(iSub,:,:,:) = switchTask3;

    % create (univariate) regressors of interest for the multivariate analyses: 
    % relational structure (-1: -Corr, 1: +Corr). 
    % pe x structure interaction.
    % Note that in both of these we're leaving the control trials with NaNs as we're not
    % interested in them - they are not part of the neural effects. 
    
    relationalStructReg     = nan(size(posCorrTrialsInds)); % we're going to leave the control trials with NaNs
    relationalStructReg(posCorrTrialsInds) = 1;
    relationalStructReg(negCorrTrialsInds) = -1;
    
    peXstructureInteractionReg = rl.STRUCT.crPe(iSub,:,:,:) .* relationalStructReg; 
    
    latentVarReg = nan(size(posCorrTrialsInds));
    latentVarReg(posCorrTrialsInds | negCorrTrialsInds) = 1;
    latentVarReg(zeroCorrTrialsInds) = -1;
    
    % do the same for vCh_STRUCT and vCh_STRUCT-vCh_NAIVE
    vChSTRUCTReg = rl.STRUCT.vCh(iSub,:,:,:);
    vChSTRUCTReg(zeroCorrTrialsInds) = nan;
    
    vChDiffReg = rl.STRUCT.vCh(iSub,:,:,:) - rl.NAIVE.vCh(iSub,:,:,:);
    vChDiffReg(zeroCorrTrialsInds) = nan;
    
    % calculate correlations with confound regressors
    stats.confoundsCorr.corrCoefs(iSub,1,1) = corr(relationalStructReg(:),subjectRts(:),'rows','complete'); % 'rows','complete' ignores NaNs
    stats.confoundsCorr.corrCoefs(iSub,1,2) = corr(relationalStructReg(:),log(subjectRts(:)),'rows','complete'); % 'rows','complete' ignores NaNs
    stats.confoundsCorr.corrCoefs(iSub,1,3) = corr(relationalStructReg(:),double(correct(:)),'rows','complete'); % 'rows','complete' ignores NaNs
    stats.confoundsCorr.corrCoefs(iSub,1,4) = corr(relationalStructReg(:),double(switchTask1(:)),'rows','complete'); % 'rows','complete' ignores NaNs
    stats.confoundsCorr.corrCoefs(iSub,1,5) = corr(relationalStructReg(:),double(switchTask2(:)),'rows','complete'); % 'rows','complete' ignores NaNs
    stats.confoundsCorr.corrCoefs(iSub,1,6) = corr(relationalStructReg(:),double(switchTask3(:)),'rows','complete'); % 'rows','complete' ignores NaNs
    stats.confoundsCorr.corrCoefs(iSub,2,1) = corr(peXstructureInteractionReg(:),subjectRts(:),'rows','complete'); % 'rows','complete' ignores NaNs
    stats.confoundsCorr.corrCoefs(iSub,2,2) = corr(peXstructureInteractionReg(:),log(subjectRts(:)),'rows','complete'); % 'rows','complete' ignores NaNs
    stats.confoundsCorr.corrCoefs(iSub,2,3) = corr(peXstructureInteractionReg(:),double(correct(:)),'rows','complete'); % 'rows','complete' ignores NaNs
    stats.confoundsCorr.corrCoefs(iSub,2,4) = corr(peXstructureInteractionReg(:),double(switchTask1(:)),'rows','complete'); % 'rows','complete' ignores NaNs
    stats.confoundsCorr.corrCoefs(iSub,2,5) = corr(peXstructureInteractionReg(:),double(switchTask2(:)),'rows','complete'); % 'rows','complete' ignores NaNs
    stats.confoundsCorr.corrCoefs(iSub,2,6) = corr(peXstructureInteractionReg(:),double(switchTask3(:)),'rows','complete'); % 'rows','complete' ignores NaNs    
    stats.confoundsCorr.corrCoefs(iSub,3,1) = corr(latentVarReg(:),subjectRts(:),'rows','complete'); % 'rows','complete' ignores NaNs
    stats.confoundsCorr.corrCoefs(iSub,3,2) = corr(latentVarReg(:),log(subjectRts(:)),'rows','complete'); % 'rows','complete' ignores NaNs
    stats.confoundsCorr.corrCoefs(iSub,3,3) = corr(latentVarReg(:),double(correct(:)),'rows','complete'); % 'rows','complete' ignores NaNs
    stats.confoundsCorr.corrCoefs(iSub,3,4) = corr(latentVarReg(:),double(switchTask1(:)),'rows','complete'); % 'rows','complete' ignores NaNs
    stats.confoundsCorr.corrCoefs(iSub,3,5) = corr(latentVarReg(:),double(switchTask2(:)),'rows','complete'); % 'rows','complete' ignores NaNs
    stats.confoundsCorr.corrCoefs(iSub,3,6) = corr(latentVarReg(:),double(switchTask3(:)),'rows','complete'); % 'rows','complete' ignores NaNs    
    stats.confoundsCorr.corrCoefs(iSub,4,1) = corr(vChSTRUCTReg(:),subjectRts(:),'rows','complete'); % 'rows','complete' ignores NaNs
    stats.confoundsCorr.corrCoefs(iSub,4,2) = corr(vChSTRUCTReg(:),log(subjectRts(:)),'rows','complete'); % 'rows','complete' ignores NaNs
    stats.confoundsCorr.corrCoefs(iSub,4,3) = corr(vChSTRUCTReg(:),double(correct(:)),'rows','complete'); % 'rows','complete' ignores NaNs
    stats.confoundsCorr.corrCoefs(iSub,4,4) = corr(vChSTRUCTReg(:),double(switchTask1(:)),'rows','complete'); % 'rows','complete' ignores NaNs
    stats.confoundsCorr.corrCoefs(iSub,4,5) = corr(vChSTRUCTReg(:),double(switchTask2(:)),'rows','complete'); % 'rows','complete' ignores NaNs
    stats.confoundsCorr.corrCoefs(iSub,4,6) = corr(vChSTRUCTReg(:),double(switchTask3(:)),'rows','complete'); % 'rows','complete' ignores NaNs    
    stats.confoundsCorr.corrCoefs(iSub,5,1) = corr(vChDiffReg(:),subjectRts(:),'rows','complete'); % 'rows','complete' ignores NaNs
    stats.confoundsCorr.corrCoefs(iSub,5,2) = corr(vChDiffReg(:),log(subjectRts(:)),'rows','complete'); % 'rows','complete' ignores NaNs
    stats.confoundsCorr.corrCoefs(iSub,5,3) = corr(vChDiffReg(:),double(correct(:)),'rows','complete'); % 'rows','complete' ignores NaNs
    stats.confoundsCorr.corrCoefs(iSub,5,4) = corr(vChDiffReg(:),switchTask1(:),'rows','complete'); % 'rows','complete' ignores NaNs
    stats.confoundsCorr.corrCoefs(iSub,5,5) = corr(vChDiffReg(:),switchTask2(:),'rows','complete'); % 'rows','complete' ignores NaNs
    stats.confoundsCorr.corrCoefs(iSub,5,6) = corr(vChDiffReg(:),switchTask3(:),'rows','complete'); % 'rows','complete' ignores NaNs    
end
% do stats on correlation coefficients
stats.confoundsCorr.P = nan(size(stats.confoundsCorr.corrCoefs,2),size(stats.confoundsCorr.corrCoefs,3));
stats.confoundsCorr.T = nan(size(stats.confoundsCorr.corrCoefs,2),size(stats.confoundsCorr.corrCoefs,3));
for iInterestReg = 1:size(stats.confoundsCorr.P,1)
    for iConfoundReg = 1:size(stats.confoundsCorr.P,2)
        [~,stats.confoundsCorr.P(iInterestReg,iConfoundReg),~,S] = ttest(atanh(squeeze(stats.confoundsCorr.corrCoefs(:,iInterestReg,iConfoundReg))));
        stats.confoundsCorr.T(iInterestReg,iConfoundReg) = S.tstat;
    end
end

%
disp('Test for differences in accuracy between different conditions') 

% +Corr vs -Corr
[~,P,~,S] = ttest(stats.accuracy.means(:,1),stats.accuracy.means(:,2),'tail','both');
stats.accuracy.tests.posVsNeg.P = P; 
stats.accuracy.tests.posVsNeg.t = S.tstat; 
fprintf('+Corr vs -Corr: P=%g, t=%g\n',P,S.tstat);

% +Corr vs 0Corr vs 
[~,P,~,S] = ttest(stats.accuracy.means(:,1),stats.accuracy.means(:,3),'tail','both');
stats.accuracy.tests.posVsZero.P = P; 
stats.accuracy.tests.posVsZero.t = S.tstat; 
fprintf('+Corr vs 0Corr: P=%g, t=%g\n',P,S.tstat);

% -Corr vs 0Corr vs 
[~,P,~,S] = ttest(stats.accuracy.means(:,2),stats.accuracy.means(:,3),'tail','both');
stats.accuracy.tests.negVsZero.P = P; 
stats.accuracy.tests.negVsZero.t = S.tstat; 
fprintf('-Corr vs 0Corr: P=%g, t=%g\n',P,S.tstat);

disp('Test for differences in log(RT) between different conditions') 

% +Corr vs -Corr
[~,P,~,S] = ttest(stats.logrt.means(:,1),stats.logrt.means(:,2),'tail','both');
stats.logrt.tests.posVsNeg.P = P; 
stats.logrt.tests.posVsNeg.t = S.tstat; 
fprintf('+Corr vs -Corr: P=%g, t=%g\n',P,S.tstat);

% +Corr vs 0Corr vs 
[~,P,~,S] = ttest(stats.logrt.means(:,1),stats.logrt.means(:,3),'tail','both');
stats.logrt.tests.posVsZero.P = P; 
stats.logrt.tests.posVsZero.t = S.tstat; 
fprintf('+Corr vs 0Corr: P=%g, t=%g\n',P,S.tstat);

% -Corr vs 0Corr vs 
[~,P,~,S] = ttest(stats.logrt.means(:,2),stats.logrt.means(:,3),'tail','both');
stats.logrt.tests.negVsZero.P = P; 
stats.logrt.tests.negVsZero.t = S.tstat; 
fprintf('-Corr vs 0Corr: P=%g, t=%g\n',P,S.tstat);


disp('Test for differences in log(RT) between different task switching groups') 

% same vs related
[~,P,~,S] = ttest(stats.logrtByTaskSwitchConds.means(:,1),stats.logrtByTaskSwitchConds.means(:,2),'tail','left');
stats.logrtByTaskSwitchCondstests.sameVsRelated.P = P; 
stats.logrtByTaskSwitchCondstests.sameVsRelated.t = S.tstat; 
fprintf('same vs related: P=%g, t=%g\n',P,S.tstat);

% same vs control vs 
[~,P,~,S] = ttest(stats.logrtByTaskSwitchConds.means(:,1),stats.logrtByTaskSwitchConds.means(:,3),'tail','left');
stats.logrtByTaskSwitchCondstests.sameVsControl.P = P; 
stats.logrtByTaskSwitchCondstests.sameVsControl.t = S.tstat; 
fprintf('same vs control: P=%g, t=%g\n',P,S.tstat);

% related vs control vs 
[~,P,~,S] = ttest(stats.logrtByTaskSwitchConds.means(:,2),stats.logrtByTaskSwitchConds.means(:,3),'tail','left');
stats.logrtByTaskSwitchCondstests.relatedVsControl.P = P; 
stats.logrtByTaskSwitchCondstests.relatedVsControl.t = S.tstat; 
fprintf('related vs control: P=%g, t=%g\n',P,S.tstat);

% # plots
% plot accuracy
figure

subplot(1,3,1)
hold on
ylabel('%correct')
h3 = boxplot(stats.accuracy.means,'symbol','r.','outliersize',14);
set(h3,'LineWidth',1)
set(gca,'YTick',[0.7:0.1:1],'FontSize',16,'XTickLabels',{'+Corr','-Corr','0Corr'},'TickLabelInterpreter','tex')

% plot logRT split by +/-/0Corr
subplot(1,3,2)
hold on
ylabel('log(rt)')
h = boxplot(stats.logrt.means,'symbol','r.','outliersize',14);
set(h,'LineWidth',1)
set(gca,'FontSize',16,'XTickLabels',{'+Corr','-Corr','0Corr'},'TickLabelInterpreter','tex')

% plot logRT  split by task switching group
subplot(1,3,3)
hold on
ylabel('log(rt)')
h = boxplot(stats.logrtByTaskSwitchConds.means,'symbol','r.','outliersize',14);
set(h,'LineWidth',1)
set(gca,'FontSize',16,'XTickLabels',{'same','related','control'},'TickLabelInterpreter','tex')

% Plot correlations of regressors of interest with confound regressors
figure;
for iRegInterest=1:5
    subplot(2,3,iRegInterest)
    hold on
    plot([0.5,length(stats.confoundsCorr.confoundRegsNames)+0.5],[0,0],'k--')
	h = boxplot(squeeze(stats.confoundsCorr.corrCoefs(:,iRegInterest,:)),'symbol','r.','outliersize',14);
    set(h,'LineWidth',1)
    set(gca,'XTickLabels',stats.confoundsCorr.confoundRegsNames,'XTick',1:length(stats.confoundsCorr.confoundRegsNames),'FontSize',16,'XTickLabelRotation',45)
    ylabel('corr coef')
    title(stats.confoundsCorr.interestRegsNames{iRegInterest})
end


save([rootData '/bhv/processed/behMatlabVars.mat'],'stats','confounds','-append');


%% Cross-validate the models cross the 4 experimental conditions
% This procedure is very similar to what we did so far, except here we
% train the data on each condition separately rather than collapse
% over blocks with the same structure, and then test data from all
% conditions. 

load([rootData '/bhv/processed/behMatlabVars.mat'],'B','B_ins','timings');

R = 20; % number of random initialisations of parameters, to ensure fitting is not local minima
theta_xVal_NAIVE  = nan(nSub,nCond,2); % 4 block types (+/-Corr, 2 stim sets). 2 parameters to fit (alpha, beta)
theta_xVal_STRUCT = nan(nSub,nCond,5); % 4 block types (+/-Corr, 2 stim sets). 5 parameters to fit (alpha, beta, 3 cross-terms)
negLogLiklihood_xVal_NAIVE = nan(nSub,nCond,nCond); % use theta from one condition to test on another condition
negLogLiklihood_xVal_STRUCT = nan(nSub,nCond,nCond); % use theta from one condition to test on another condition

for iSub = 1:length(subjects)
    sub = subjects{iSub};
    % we will be concatenating data from the single pre-scan block and the 2 scanning
    % blocks (from the 2 runs) -  initialise this now.
    cues_cat     = nan(nCond,nTrialsPerBlock_instruct + 2*nTrialsPerBlock);
    outcomes_cat = nan(nCond,nTrialsPerBlock_instruct + 2*nTrialsPerBlock);
    choices_cat  = nan(nCond,nTrialsPerBlock_instruct + 2*nTrialsPerBlock);
    
    % Get a vector of trial numbers of last trial of blocks - model
    % estimations should reset at the beginning of blocks. should be
    % [42, 72, 102]
    blocks = [nTrialsPerBlock_instruct,nTrialsPerBlock_instruct + nTrialsPerBlock, nTrialsPerBlock_instruct + 2*nTrialsPerBlock];
    
    for iCond=1:nCond; % loop over conditions
        % concatenating data from the single pre-scan block and the 2 scanning blocks (from
        % the 2 runs)
        cues_cat(iCond,:)     = [squeeze(B_ins.cues(iSub,1,iCond,:))',squeeze(B.cues(iSub,1,iCond,:))',squeeze(B.cues(iSub,2,iCond,:))'];
        outcomes_cat(iCond,:) = [squeeze(B_ins.outcomes(iSub,1,iCond,:))',squeeze(B.outcomes(iSub,1,iCond,:))',squeeze(B.outcomes(iSub,2,iCond,:))'];
        choices_cat(iCond,:)  = [squeeze(B_ins.choices(iSub,1,iCond,:))',squeeze(B.choices(iSub,1,iCond,:))',squeeze(B.choices(iSub,2,iCond,:))'];
        
        % store maximum likelihood parameters for NAIVE: [alpha,beta].
        [theta_xVal_NAIVE(iSub,iCond,:), ~] = fitRlParams(cues_cat(iCond,:),outcomes_cat(iCond,:),choices_cat(iCond,:),false,blocks,R);
        % store maximum likelihood parameters for STRUCT: [alpha,beta,H_AB,H_AC,H_BC].
        [theta_xVal_STRUCT(iSub,iCond,:), ~] = fitRlParams(cues_cat(iCond,:),outcomes_cat(iCond,:),choices_cat(iCond,:),true,blocks,R);
    end
    
    % get cross validated log likelihoods
    for iCond=1:nCond % condition of test data
        for jCond=1:nCond % condition of training data
            [negLogLiklihood_xVal_NAIVE(iSub,iCond,jCond),~,~]  = calcRlLogLikelihood(theta_xVal_NAIVE(iSub,jCond,:),cues_cat(iCond,:),outcomes_cat(iCond,:),choices_cat(iCond,:),false,blocks);
            [negLogLiklihood_xVal_STRUCT(iSub,iCond,jCond),~,~] = calcRlLogLikelihood(theta_xVal_STRUCT(iSub,jCond,:),cues_cat(iCond,:),outcomes_cat(iCond,:),choices_cat(iCond,:),true,blocks);
        end
    end
end
xVal.thetaMl.NAIVE    = theta_xVal_NAIVE;
xVal.thetaMl.STRUCT   = theta_xVal_STRUCT;
xVal.negLogLik.NAIVE  = negLogLiklihood_xVal_NAIVE;
xVal.negLogLik.STRUCT = negLogLiklihood_xVal_STRUCT;
save([rootData '/bhv/processed/behMatlabVars.mat'],'xVal','-append');


%% cross-validate trial-by-trial values of RL variables (STRUCT and NAIVE models)
% using parameters fitted on different conditions

% RL variables that are tracked (for both NAIVE and STRUCT models, and for
% each of the nCond x nCond combinations of training and testing data):

% g      : Outcome probability estimation.
% vCh    : Value of the chosen action, which always has the same magnitude
%          but has an opposite sign to g on rejected trials
% epsilon: Outcome prediction error
% crPe   : The "correctness" prediction error, which has the same magnitude
%           as epsilon but a sign that depends on the congruence of the
%           outcome and the subject's choice.
load([rootData '/bhv/processed/behMatlabVars.mat'],'xVal','B','B_ins','rl')

xVal.rl.NAIVE.g       = nan(nSub,nRun,nCond,nCond,nTrialsPerBlock);
xVal.rl.NAIVE.vCh     = nan(nSub,nRun,nCond,nCond,nTrialsPerBlock);
xVal.rl.NAIVE.epsilon = nan(nSub,nRun,nCond,nCond,nTrialsPerBlock);
xVal.rl.NAIVE.crPe    = nan(nSub,nRun,nCond,nCond,nTrialsPerBlock);

xVal.rl.STRUCT.g       = nan(nSub,nRun,nCond,nCond,nTrialsPerBlock);
xVal.rl.STRUCT.vCh     = nan(nSub,nRun,nCond,nCond,nTrialsPerBlock);
xVal.rl.STRUCT.epsilon = nan(nSub,nRun,nCond,nCond,nTrialsPerBlock);
xVal.rl.STRUCT.crPe    = nan(nSub,nRun,nCond,nCond,nTrialsPerBlock);

% cues are coded as 0,1,2; for Matlab indexing move to be 1,2,3
cues = B.cues+1;

for iSub=1:nSub
    for iRun=1:nRun
        for iCond=1:nCond % testing data
            for jCond=1:nCond % training data
                
                % get parameters from training data - jCond
                alpha_NAIVE  = squeeze(xVal.thetaMl.NAIVE(iSub,jCond,1));
                alpha_STRUCT = squeeze(xVal.thetaMl.STRUCT(iSub,jCond,1));
                crossTerms   = squeeze(xVal.thetaMl.STRUCT(iSub,jCond,3:5))';
                
                % initialise g at 0 for all three stimuli
                g_NAIVE = zeros(3,1); % 3 stimuli
                g_STRUCT = zeros(3,1); % 3 stimuli
                
                for t=1:nTrialsPerBlock
                    % find current cue from test data - iCond. Cues are coded as 0,1,2, so need to add
                    % +1 for Matlab indexing. 
                    c = B.cues(iSub,iRun,iCond,t) +1;
                    
                    % calculate outcome probability etimation and outcome
                    % prediciton error
                    xVal.rl.NAIVE.g(iSub,iRun,iCond,jCond,t)        = g_NAIVE(c);
                    xVal.rl.NAIVE.epsilon(iSub,iRun,iCond,jCond,t)  = B.outcomes(iSub,iRun,iCond,t) - g_NAIVE(c);
                    xVal.rl.STRUCT.g(iSub,iRun,iCond,jCond,t)       = g_STRUCT(c);
                    xVal.rl.STRUCT.epsilon(iSub,iRun,iCond,jCond,t) = B.outcomes(iSub,iRun,iCond,t) - g_STRUCT(c);
                    
                    % flip sign if subject rejected for vCh and crPe.
                    if B.choices(iSub,iRun,iCond,t)==1 % trial accepted
                        xVal.rl.NAIVE.vCh(iSub,iRun,iCond,jCond,t)   = xVal.rl.NAIVE.g(iSub,iRun,iCond,t);
                        xVal.rl.NAIVE.crPe(iSub,iRun,iCond,jCond,t)  = xVal.rl.NAIVE.epsilon(iSub,iRun,iCond,t);
                        xVal.rl.STRUCT.vCh(iSub,iRun,iCond,jCond,t)  = xVal.rl.STRUCT.g(iSub,iRun,iCond,t);
                        xVal.rl.STRUCT.crPe(iSub,iRun,iCond,jCond,t) = xVal.rl.STRUCT.epsilon(iSub,iRun,iCond,t);
                    elseif B.choices(iSub,iRun,iCond,t)==0 % trial rejected
                        xVal.rl.NAIVE.vCh(iSub,iRun,iCond,jCond,t)   = - xVal.rl.NAIVE.g(iSub,iRun,iCond,t);
                        xVal.rl.NAIVE.crPe(iSub,iRun,iCond,jCond,t)  = - xVal.rl.NAIVE.epsilon(iSub,iRun,iCond,t);
                        xVal.rl.STRUCT.vCh(iSub,iRun,iCond,jCond,t)  = - xVal.rl.STRUCT.g(iSub,iRun,iCond,t);
                        xVal.rl.STRUCT.crPe(iSub,iRun,iCond,jCond,t) = - xVal.rl.STRUCT.epsilon(iSub,iRun,iCond,t);
                    end
                    g_NAIVE  = getRlValue(g_NAIVE,B.outcomes(iSub,iRun,iCond,t),B.cues(iSub,iRun,iCond,t),3,alpha_NAIVE,[0,0,0]);
                    g_STRUCT = getRlValue(g_STRUCT,B.outcomes(iSub,iRun,iCond,t),B.cues(iSub,iRun,iCond,t),3,alpha_STRUCT,crossTerms);
                end
            end
        end
    end
end

% obtain xVal.rl.STRUCT.g_gnrl: same dims as rl.STRUCT.g only in each
% condition using the parameters from the model trained on the different
% stimuli set, same structure condition. 

xVal.rl.STRUCT.g_gnrl = nan(size(rl.STRUCT.g));
xVal.rl.STRUCT.g_gnrl(:,:,1,:) = xVal.rl.STRUCT.g(:,:,1,2,:); % data from cond 1, training on cond 2
xVal.rl.STRUCT.g_gnrl(:,:,2,:) = xVal.rl.STRUCT.g(:,:,2,1,:); % data from cond 2, training on cond 1
xVal.rl.STRUCT.g_gnrl(:,:,3,:) = xVal.rl.STRUCT.g(:,:,3,4,:); % data from cond 3, training on cond 4
xVal.rl.STRUCT.g_gnrl(:,:,4,:) = xVal.rl.STRUCT.g(:,:,4,3,:); % data from cond 4, training on cond 3

% save
save([rootData '/bhv/processed/behMatlabVars.mat'],'xVal','-append');


%% stats
% training and testing on same data
LL_diagSTRUCT = mean([xVal.negLogLik.STRUCT(:,1,1),xVal.negLogLik.STRUCT(:,2,2),xVal.negLogLik.STRUCT(:,3,3),xVal.negLogLik.STRUCT(:,4,4)],2);
LL_diagNAIVE  = mean([xVal.negLogLik.NAIVE(:,1,1),xVal.negLogLik.NAIVE(:,2,2),xVal.negLogLik.NAIVE(:,3,3),xVal.negLogLik.NAIVE(:,4,4)],2);

% Cross validated - within and between structures (conds 1,2: -Corr; conds
% 3,4: +Corr)

LL_xValSameStruct = mean([xVal.negLogLik.STRUCT(:,1,2),xVal.negLogLik.STRUCT(:,2,1),xVal.negLogLik.STRUCT(:,3,4),xVal.negLogLik.STRUCT(:,4,3)],2);
LL_xValDiffStruct = mean([xVal.negLogLik.STRUCT(:,1,3),xVal.negLogLik.STRUCT(:,1,4),xVal.negLogLik.STRUCT(:,2,3),xVal.negLogLik.STRUCT(:,2,4),...
    xVal.negLogLik.STRUCT(:,3,1),xVal.negLogLik.STRUCT(:,4,1),xVal.negLogLik.STRUCT(:,4,1),xVal.negLogLik.STRUCT(:,4,2)],2);

% difference between the within-structure x-val STRUCT model and STRUCT model
% trained and tested on the same data.

[~,P,~,S] = ttest(LL_diagSTRUCT,LL_xValSameStruct);
fprintf('within-struct x-Val vs STRUCT diag: P=%g, t=%g\n',P,S.tstat);

% difference between the within-structure x-val STRUCT model and NAIVE model
% trained and tested on the same data.

[~,P,~,S] = ttest(LL_diagNAIVE,LL_xValSameStruct);
fprintf('within-struct x-Val vs NAIVE diag: P=%g, t=%g\n',P,S.tstat);

% difference between the within-structure x-val STRUCT model and
% between-structures x-val STRUCT model.

[~,P,~,S] = ttest(LL_xValSameStruct,LL_xValDiffStruct);
fprintf('within-struct x-Val vs between-struct x-Val: P=%g, t=%g\n',P,S.tstat);


%% plots - cross-validaeted models (figs 2a and 2b).

% cross-validated log-likelihood matrices: 
% The within-structure cross-validated STRUCT model explains the data
% better then the NAIVE model, even when the latter is tested on the same
% data it was trained on. 

figure; 
% sum log likelihood across subjects
ml_NAIVE = squeeze(sum(xVal.negLogLik.NAIVE,1)); 
ml_STRUCT = squeeze(sum(xVal.negLogLik.STRUCT,1));
% plot both matrices in the same scale.
colormap hot
condStr = {'-1','-2','+1','+2'};

subplot(1,2,1)
imagesc(ml_NAIVE,[min([min(min(ml_NAIVE)),min(min(ml_STRUCT))]),max([max(max(ml_NAIVE)),max(max(ml_STRUCT))])])
set(gca,'YTick',1:4,'YTickLabel',condStr,'FontSize',16,'XTick',1:4,'XTickLabel',condStr)
xlabel('test'); ylabel('train')
title('NAIVE')
axis square

subplot(1,2,2)
imagesc(ml_STRUCT,[min([min(min(ml_NAIVE)),min(min(ml_STRUCT))]),max([max(max(ml_NAIVE)),max(max(ml_STRUCT))])])
set(gca,'YTick',1:4,'YTickLabel',condStr,'FontSize',16,'XTick',1:4,'XTickLabel',condStr)
xlabel('test'); ylabel('train')
title('STRUCT')
axis square

% ###
% plot histograms of g for accepted and rejecteds trials, particularly in
% trials where STRUCT and NAIVe models disagree. 

% find indeces of trials where STRUCT and NAIVE make different choice
% predictions

I = sign(xVal.rl.STRUCT.g_gnrl(:))~=sign(rl.NAIVE.g(:));

figure

subplot(1,2,1)
hold on
title('diff trials')
tmp1 = xVal.rl.STRUCT.g_gnrl(:);
tmp2 = B.choices(:);
histogram(tmp1(tmp2==1 & I))
histogram(tmp1(tmp2==0 & I))
legend('accept','reject')
set(gca,'FontSize',20,'YTick',[])
xlabel('g STRUCT x-Val')

subplot(1,2,2)
hold on
title('diff trials')
tmp1 = rl.NAIVE.g(:);
tmp2 = B.choices(:);
histogram(tmp1(tmp2==1 & I))
histogram(tmp1(tmp2==0 & I))
legend('accept','reject')
set(gca,'FontSize',20,'YTick',[])
xlabel('g NAIVE')

%%


