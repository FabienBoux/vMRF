
addpath(genpath(fullfile(pwd, 'functions')))

%% Regression approach Vs Grid search approach

% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dico_regression     = 'files/16-07-2018_dico_rand_4params_160000signals.mat';
dico_regression     = 'files/dico_rand_Ngrand_160000signals.mat';

signal_test_number  = 10000;
cstr.Sigma          = '';
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load(dico_regression)

% For each parameter size
parfor iter = 1:size(Y,2)

    fprintf('%d \n', iter)

    % Build test dataset
    rand_perm           = randperm(length(X), signal_test_number);
    Xtest{iter}         = abs(X(rand_perm, :));
    Ytest{iter}         = Y(rand_perm, 1:iter);

    % Build learning dataset
    vect                = ones(size(X,1),1); vect(rand_perm) = 0;
    Xtrain{iter}        = X(logical(vect), :);
    Ytrain{iter}        = Y(logical(vect), 1:iter);
    dico_size(iter)     = size(Ytrain{iter},1);

    % Train model
    [theta{iter}, ~]  	= EstimateInverseFunction(Ytrain{iter}, Xtrain{iter}, 20, 5-iter, 100, cstr, 0);

    % Predict using the trained model
    Ypredict{iter}      = EstimateParametersFromModel(Xtest{iter}, theta{iter}); 
%     for i = 1:iter
%         Ypredict{iter}(Ypredict{iter}(:,i) > max(Y(:,i))) = nan;
%         Ypredict{iter}(Ypredict{iter}(:,i) < min(Y(:,i))) = nan;
%     end

    % Compute estimation error
    [Rmse{iter}, Nrmse{iter}, Mae{iter}] = EvaluateEstimation(Ytest{iter}, Ypredict{iter}, Ytrain{iter});
end




