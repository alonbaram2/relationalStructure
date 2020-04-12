clear
rootScripts = '/home/fs0/abaram/scripts/NN2020';
rootData    = '/home/fs0/abaram/scratch/NN2020_noRound';
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
B.correct  = nan(nSub,nRun,nCond,nTrialsPerBlock); 
B.rt       = nan(nSub,nRun,nCond,nTrialsPerBlock); 

timings.stim    = nan(nSub,nRun,nCond,nTrialsPerBlock);
timings.outcome = nan(nSub,nRun,nCond,nTrialsPerBlock);
timings.choice  = nan(nSub,nRun,nCond,nTrialsPerBlock);

B_ins.cues     = nan(nSub,nRun,nCond,nTrialsPerBlock_instruct); 
B_ins.outcomes = nan(nSub,nRun,nCond,nTrialsPerBlock_instruct); 
B_ins.choices  = nan(nSub,nRun,nCond,nTrialsPerBlock_instruct);
B_ins.correct  = nan(nSub,nRun,nCond,nTrialsPerBlock_instruct); 
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
            B.correct(iSub,run,cond,:)  = data.run{run}.cond{cond}.correct; % congruency between choice and outcome.            
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
        B_ins.correct(iSub,run,cond,:)  = data.run{run}.cond{cond}.correct; % congruency between choice and outcome.
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
fprintf('Test cross-terms Corr vs 0Corr: P=%g, t=%g',P,S.tstat);

% ###
% model comparison
N = negLogLik.N; % number of examples (trials)

disp('model comparison: smaller is better')
% NAIVE
LL = sum(sum( - negLogLik.NAIVE)); % LL: log likelihood. collapse over -/+ Corr blocks and over subjects
k = 2 * nSub; %number of parameters fitted (2 per subject: alpha and beta)
AIC = 2*k - 2*LL;
BIC = log(N)*k - 2*LL;
disp(sprintf('NAIVE: -LL: %g , AIC: %g ,  BIC: %g' ,-LL,AIC,BIC))
stats.modelComparison.NAIVE.negLL = -LL;
stats.modelComparison.NAIVE.AIC = AIC;
stats.modelComparison.NAIVE.BIC = BIC;

% STRUCT
LL = sum(sum( - negLogLik.STRUCT)); % LL: log likelihood. collapse over -/+ Corr blocks and over subjects
k = 5 * nSub; %number of parameters fitted (5 per subject: alpha, beta and 3 cross terms)
AIC = 2*k - 2*LL;
BIC = log(N)*k - 2*LL;
disp(sprintf('STRUCT: -LL: %g , AIC: %g ,  BIC: %g' ,-LL,AIC,BIC))
stats.modelComparison.STRUCT.negLL = -LL;
stats.modelComparison.STRUCT.AIC = AIC;
stats.modelComparison.STRUCT.BIC = BIC;

% ###
% plot histograms of g for accepted and rejecteds trials, particularly in
% trials where STRUCT and NAIVe models disagree. 

% find indeces of trials where STRUCT and NAIVE make different choice
% predictions

I = sign(rl.STRUCT.g(:))~=sign(rl.NAIVE.g(:));

figure
subplot(2,2,1)
hold on
title('all trials')
histogram(rl.STRUCT.g(B.choices==1))
histogram(rl.STRUCT.g(B.choices==0))
legend('accept','reject')
set(gca,'FontSize',20,'YTick',[])
xlabel('g STRUCT')

subplot(2,2,2)
hold on
title('all trials')
histogram(rl.NAIVE.g(B.choices==1))
histogram(rl.NAIVE.g(B.choices==0))
legend('accept','reject')
set(gca,'FontSize',20,'YTick',[])
xlabel('g NAIVE')

subplot(2,2,3)
hold on
title('diff trials')
tmp1 = rl.STRUCT.g(:);
tmp2 = B.choices(:);
histogram(tmp1(tmp2==1 & I))
histogram(tmp1(tmp2==0 & I))
legend('accept','reject')
set(gca,'FontSize',20,'YTick',[])
xlabel('g STRUCT')

subplot(2,2,4)
hold on
title('diff trials')
tmp1 = rl.NAIVE.g(:);
tmp2 = B.choices(:);
histogram(tmp1(tmp2==1 & I))
histogram(tmp1(tmp2==0 & I))
legend('accept','reject')
set(gca,'FontSize',20,'YTick',[])
xlabel('g NAIVE')

% #############
%
% Accuracy and reaction time analyses
%

load([rootData '/bhv/processed/behMatlabVars.mat'],'B','rl')


% Accuracy (i.e. choosing correctly according to STRUC model): Is there a
% difference in between different conditions?

stats.accuracy.means = nan (nSub,3); % percent correct for 3 conditions: [+,-,0] Corr
stats.logrt.means = nan (nSub,3); % mean log(rt) for 3 conditions: [+,-,0] Corr

for iSub=1:length(subjects)
    % find the correct choice trials according to STRUCT model    
    correct_STRUCT =  sign(B.choices(iSub,:,:,:).*2-1) == sign(rl.STRUCT.g(iSub,:,:,:)); % change choices to be {-1,1} instead of {0,1}
    % store rt of specific subject in a temporary variable (done just for
    % simplicity of coding)
    tmp_subjectRts = B.rt(iSub,:,:,:);
        
    
    % find indeces of different trial types;  %
    % find indeces of all +Corr blocks:    
    posCorrBlocksInds = false(size(B.cues(iSub,:,:,:)));
    posCorrBlocksInds(1,:,3:4,:) = true;
    % find all A and B cues in +Corr blocks
    posCorrTrialsInds = posCorrBlocksInds & (B.cues(iSub,:,:,:) ~= 2 );  
    % percent correct in +Corr trials
    stats.accuracy.means(iSub,1) = sum(correct_STRUCT(posCorrTrialsInds)) ./ length(correct_STRUCT(posCorrTrialsInds));
    % find mean log(rt) in +Corr
    stats.logrt.means(iSub,1) = mean(tmp_subjectRts(posCorrTrialsInds));
       
    % First find indeces of all -Corr blocks:    
    negCorrBlocksInds = false(size(B.cues(iSub,:,:,:)));
    negCorrBlocksInds(1,:,1:2,:) = true;
    % find all A and B cues in -Corr blocks
    negCorrTrialsInds = negCorrBlocksInds & (B.cues(iSub,:,:,:) ~= 2 );
    % percent correct in -Corr trials
    stats.accuracy.means(iSub,2) = sum(correct_STRUCT(negCorrTrialsInds)) ./ length(correct_STRUCT(negCorrTrialsInds));
    % find mean log(rt) in -Corr
    stats.logrt.means(iSub,2) = mean(tmp_subjectRts(negCorrTrialsInds));
    
    % find all C cues
    zeroCorrTrialsInds = false(size(B.cues(iSub,:,:,:)));
    zeroCorrTrialsInds(B.cues(iSub,:,:,:)==2) = true;
    % percent correct in 0Corr trials
    stats.accuracy.means(iSub,3) = sum(correct_STRUCT(zeroCorrTrialsInds)) ./ length(correct_STRUCT(zeroCorrTrialsInds));
    % find mean log(rt) in 0Corr
    stats.logrt.means(iSub,3) = mean(tmp_subjectRts(zeroCorrTrialsInds));
end

%
disp('Test for differences in accuracy between different conditions') 

% +Corr vs -Corr
[~,P,~,S] = ttest(stats.accuracy.means(:,1),stats.accuracy.means(:,2),'tail','both');
stats.accuracy.tests.posVsNeg.P = P; 
stats.accuracy.tests.posVsNeg.t = S.tstat; 
fprintf('+Corr vs -Corr: P=%g, t=%g',P,S.tstat);

% +Corr vs 0Corr vs 
[~,P,~,S] = ttest(stats.accuracy.means(:,1),stats.accuracy.means(:,3),'tail','both');
stats.accuracy.tests.posVsZero.P = P; 
stats.accuracy.tests.posVsZero.t = S.tstat; 
fprintf('+Corr vs 0Corr: P=%g, t=%g',P,S.tstat);

% -Corr vs 0Corr vs 
[~,P,~,S] = ttest(stats.accuracy.means(:,2),stats.accuracy.means(:,3),'tail','both');
stats.accuracy.tests.negVsZero.P = P; 
stats.accuracy.tests.negVsZero.t = S.tstat; 
fprintf('-Corr vs 0Corr: P=%g, t=%g',P,S.tstat);

% plot accuracy
figure
hold on
title('%correct')
h3 = boxplot(stats.accuracy.means);
set(h3,'LineWidth',2)
set(gca,'YTick',[0.7:0.1:1],'FontSize',16,'XTickLabels',{'+Corr','-Corr','0Corr'},'TickLabelInterpreter','tex')


disp('Test for differences in log(RT) between different conditions') 

% +Corr vs -Corr
[~,P,~,S] = ttest(stats.logrt.means(:,1),stats.logrt.means(:,2),'tail','both');
stats.logrt.tests.posVsNeg.P = P; 
stats.logrt.tests.posVsNeg.t = S.tstat; 
fprintf('+Corr vs -Corr: P=%g, t=%g',P,S.tstat);

% +Corr vs 0Corr vs 
[~,P,~,S] = ttest(stats.logrt.means(:,1),stats.logrt.means(:,3),'tail','both');
stats.logrt.tests.posVsZero.P = P; 
stats.logrt.tests.posVsZero.t = S.tstat; 
fprintf('+Corr vs 0Corr: P=%g, t=%g',P,S.tstat);

% -Corr vs 0Corr vs 
[~,P,~,S] = ttest(stats.logrt.means(:,2),stats.logrt.means(:,3),'tail','both');
stats.logrt.tests.negVsZero.P = P; 
stats.logrt.tests.negVsZero.t = S.tstat; 
fprintf('-Corr vs 0Corr: P=%g, t=%g',P,S.tstat);

% plot logRT
figure
hold on
title('log(rt)')
h4 = boxplot(stats.logrt.means);
set(h4,'LineWidth',2)
set(gca,'FontSize',16,'XTickLabels',{'+Corr','-Corr','0Corr'},'TickLabelInterpreter','tex')

save([rootData '/bhv/processed/behMatlabVars.mat'],'stats','-append');