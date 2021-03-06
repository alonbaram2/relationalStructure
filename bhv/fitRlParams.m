function [theta_ML, neglogLik_ML] = fitRlParams(cues,outcomes,choices,STRUCT_flag,blocks,R)
% cues: Tx1 vecto of stimuli (1: stim A; 2: stim B; 3: stim C)
% outcomes: Tx1 vector of outcomes (0: bad outcome; 1: good outcome. 
% choices: Tx1 vector of participant/agent choices: 0: reject, 1: accept
% STRUCT_flag: 1: STRUCT model; 0: NAIVE model. i.e. fit cross-terms?

NeglogLik_ML_all = nan(R,1);
if STRUCT_flag
    theta_ML_all= nan(R,5); % [alpha, beta, C_ab, C_ac, C_bc]
else
    theta_ML_all= nan(R,2); % [alpha, beta]
    
end

LB = [0,0]; % Lower bound for [alpha, beta]
UB = [1,8]; % Upper bound for [alpha, beta]
X0 = [LB(1) + (UB(1)-LB(1)).*rand(R,1) , LB(2) + (UB(2)-LB(2)).*rand(R,1)]; % random initialisation

if STRUCT_flag
    LB = [LB, -1, -1, -1]; % -1 is lower bound for cross-terms
    UB = [UB, 1, 1, 1]; % 1 is upper  bound for cross-terms
    X0 = [X0, -1 + 2*rand(R,1),-1 + 2*rand(R,1),-1 + 2*rand(R,1)]; %  random initialisation for cross-terms
end

energy = @calcRlLogLikelihood;
OPTIM_options = optimset('Display', 'off') ;

for i = 1:R
    [theta_ML_all(i,:), NeglogLik_ML_all(i)] = fmincon(@(theta) energy(theta,cues,outcomes,choices,STRUCT_flag,blocks),X0(i,:)',[],[],[],[],LB', UB',[],OPTIM_options);
end
[neglogLik_ML,I] = min(NeglogLik_ML_all);
theta_ML = theta_ML_all(I,:);

