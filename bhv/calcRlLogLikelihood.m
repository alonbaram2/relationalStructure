function [negLogLik,VV,PP] = calcRlLogLikelihood(theta,cues,outcomes,choices,STRUCT_flag,blocks)
% throughout this script the word "value" is used to describe the model
% estimation of the probability of a "good" outcome, Note that this
% "probability" is betweeen -1 and 1, and is turned into actual probability
% (between 0 and 1) by passing it through a sigmoid function. The reason
% for the this transformation is the need for the indifference point to be
% 0 rather than 0.5 for the equation governing the use of the cross-terms
% to work properly. 

nCues = length(unique(cues)); % number of possible cues: 3
alpha= theta(1); % learning rate for experienced rewards
beta = theta(2);
if STRUCT_flag 
    C_ab = theta(3);
    C_ac = theta(4);
    C_bc = theta(5);
else % NAIVE model
    C_ab = 0;
    C_ac = 0;
    C_bc = 0;
end
corr_params = [C_ab, C_ac, C_bc];

v0 = 0; % initial value at beginning of block
v = nan(nCues,1);
v(:)=v0;

T = length(cues); % number of trials
PP = nan(T,1); % track prob of accepting through time 
                % (pass value after last trial through sigmoid).  this is
                % between 0 and 1). This is used before viewing the current
                % outcome. 
VV = nan(T,3); % track value through time,this is between -1 and 1. This is the updated value, after viewing the outcome.

for t=1:T % trials    
    if (cues(t)==0)
        c = 1;
    elseif (cues(t)==1)
        c = 2;
    elseif (cues(t)==2)
        c = 3;
    end
    
    % calculate probability of accepting in current trial, by passing the
    % updated value from the previous trial through a sigmoid) 
    p   = 1/(1+exp(- beta*(v(c)))); % pass through sigmoid to obtain probabilities between 0 and 1. 
    PP(t) = p;
    
    % update value according to current trial outcome
    v = getRlValue(v,outcomes(t),cues(t),nCues,alpha,corr_params);
    
    % reset value before beginnings of blocks
    if any(t==blocks) 
        v(:)=v0;
    end
    VV(t,:) = v;     
end
negLogLik = sum(-log(1-PP(choices==0)))+sum(-log(PP(choices==1)));

