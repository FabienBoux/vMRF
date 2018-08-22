
load('dictionaries/multidico/2018-07-25-09:35_32-samples_5-parameters_80000-signals.mat')

l = ceil(2 * size(X,1) /3);

Xtrain  = X(1:l,:);
Ytrain  = Y(1:l,:);
Xtest   = X(l+1:end,:);
Ytest   = Y(l+1:end,:);

Kmax    = 10;

[Kbest, Kopti] = FindOptimalK(Xtrain, Ytrain, Xtest, Ytest, Kmax);