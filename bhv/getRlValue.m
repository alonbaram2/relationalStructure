function v = getRlValue(v,outcome,cue,nCues,alpha,corr_params)
C_ab = corr_params(1);
C_ac = corr_params(2);
C_bc = corr_params(3);

% find current cue. Cues are coded as 0,1,2, so need to add
% +1 for Matlab indexing
for i = 0:nCues-1
    if (cue==i)
        c = i+1;
    end
end

dv = outcome - v(c);
v(c) = v(c) + alpha*dv;

if c==1
    v(2) = (1 - abs(C_ab)*alpha) * v(2) + C_ab*alpha*outcome;
    v(3) = (1 - abs(C_ac)*alpha) * v(3) + C_ac*alpha*outcome;
elseif c==2
    v(1) = (1 - abs(C_ab)*alpha) * v(1) + C_ab*alpha*outcome;
    v(3) = (1 - abs(C_bc)*alpha) * v(3) + C_bc*alpha*outcome;
elseif c==3
    v(1) = (1 - abs(C_ac)*alpha) * v(1) + C_ac*alpha*outcome;
    v(2) = (1 - abs(C_bc)*alpha) * v(2) + C_bc*alpha*outcome;
end
    



