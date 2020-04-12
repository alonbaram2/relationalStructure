function [S, L] = mvpa_svm_OneVsAll(X, y)

% Load data stored in arrays X, y
%load('X.mat'); % examples by features matrix
%load('y.mat'); % examples by class vector 

% Setup the parameters
num_labels = 8;                % number of categories to which to classify
m = size(X, 1);               % number of examples

%fprintf('\nTraining One-vs-All Support vector machine...\n')

% loop around each example to train without it and then predict it

S=NaN(m, num_labels);
L=NaN(m, num_labels);

for i = 1:m    
    Xtrain = X; ytrain = y;
    Xtrain(i,:)=[]; ytrain(i)=[];
    for c = 1:num_labels
        
        svm_model = fitcsvm(Xtrain,(ytrain==c));
        [label,score] = predict(svm_model,X(i,:));
        S(i,c)=score(2);
        L(i,c)=label;
    end
end;

end