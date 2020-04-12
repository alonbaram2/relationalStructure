function out=calcRDMs(Y,SPM,condition)
% Y is usually full(X(T.LI{i},:)') from rsa.runSearchlight.m: raw data,
% only withing the current searchlight.

% spatial prewhitening: 
% Do GLM and Get prewhitened beta weights
B=real(rsa.spm.noiseNormalizeBeta(Y,SPM));
%take only betas of relevant conditions
B=B(logical(condition),:);
% calculate distance
RDM = pdist(B,'correlation');    
% meanDist=mean(RDM,2);           % Calulate the mean distance
RDM=RDM';
% out = [RDM(:);meanDist(:)];     % Arrange the outputs in a vector, as they are written to files
out = RDM(:);     % Arrange the outputs in a vector, as they are written to files
