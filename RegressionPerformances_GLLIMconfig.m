
addpath(genpath(fullfile(pwd, 'functions')))


%% Better K ? K = 20
% bic is not working

% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_dico = 'files/dico_random_92559_signals.mat';
signal_test_number = 1000;
signal_train_number = 5000;
snr_values = inf;
K = 1:50;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(file_dico)

Rmse = zeros(length(K), 2); Nrmse = Rmse;


for iter = 1:length(K)

    fprintf('K value %d\n', K(iter))

    rand_perm   = randperm(length(X), signal_train_number);

    Xtrain      = abs(X(rand_perm, :));
    Ytrain      = Y(rand_perm, [1 3]);

    Xtest       = abs(X(1:signal_test_number, :));
    Ytest       = Y(1:signal_test_number, [1 3]);

    Xtrain      = AddNoise(Xtrain, snr_values, 0);
    Xtest       = AddNoise(Xtest , snr_values, 0);

    [theta, r]  = EstimateInverseFunction(Ytrain, Xtrain, K(iter), 0, 30, 0);
    Ypredict  	= EstimateParametersFromModel(Xtest, theta);

    [Rmse(iter, :), Nrmse(iter, :)] = EvaluateEstimation(Ytest, Ypredict, Ytrain);
    
    D           = size(Ytrain, 2);
    L           = size(X, 2);
    N           = size(Xtrain, 1);
    bic(iter)   = -2 * sum( log(sum(r,1)) ) + K(iter) * (D*(L+1) + 0.5*L*(L+3) + 1) * log(N);
end

%% Plotting

figure
subplot(1,2,1); plot(K, mean(Nrmse,2)', 'o-', 'linewidth',1.5)
subplot(1,2,2); plot(K, bic, 'o-', 'linewidth',1.5)
    
    
    
    
    