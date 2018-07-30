
addpath(genpath(fullfile(pwd, 'functions')))

%% Regression approach Vs Grid search approach
%  Nrmse = 1 means that we estimate by the mean value


% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_dico  = 'files/dico_random_220900signals';

dico_size = [1000 2000 5000 10000:10000:100000 120000:20000:220000];
signal_test_number = 10000;

repetition = 3;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(file_dico)

for rep = 1:repetition
    for iter = 1:length(dico_size)

        rand_perm 	= randperm(length(X), dico_size(iter));
        Xtrain      = abs(X(rand_perm, :));
        Ytrain  	= Y(rand_perm, :);

        rand_perm 	= randperm(length(X), signal_test_number);
        Xtest   	= abs(X(rand_perm, :));
        Ytest     	= Y(rand_perm, :);

        Ypredict    = mean(Ytest) .* ones(size(Ytest));

        [Rmse(iter, :, rep), Nrmse(iter, :, rep)] = EvaluateEstimation(Ytest, Ypredict, Ytrain);
    end
end

%% Plotting

figure
semilogx(dico_size, mean(Nrmse,3), 'o-', 'linewidth', 1.5)
xlabel('Dictionary size'); ylabel('NRMSE');











