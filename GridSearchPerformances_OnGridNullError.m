
addpath(genpath(fullfile(pwd, 'functions')))


%% Verify that signals on grid are correctly estimate
% The error sum is null = it's ok = every signal on the grid is estimate by
% itself

% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_dico = 'files/dico_grid_216931signals.mat';

dico_size = [2000:2000:10000 15000 20000:10000:200000];
signal_test_number = 1000;

repetition = 1;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(file_dico)

% For each iteration
for rep = 1:repetition
    
    % For each dictionary size
    for iter = 1:length(dico_size)

        fprintf('Dictionary %d / %d \t %d signals \n', iter, length(dico_size), dico_size(iter))

        rand_perm   = randperm(length(X), dico_size(iter));

        Xtrain      = abs(X(rand_perm, :));
        Ytrain      = Y(rand_perm, :);

        Xtest       = abs(Xtrain(1:signal_test_number, :));
        Ytest       = Ytrain(1:signal_test_number, :);
        
        Ypredict  	= EstimateParametersFromGrid(Xtest, Xtrain, Ytrain);

        [Rmse(iter, :, rep), Nrmse(iter, :, rep)] = EvaluateEstimation(Ytest, Ypredict, Ytrain);

    end
end

%% Display results

fprintf('Sum of errors: %f \n', sum(sum(sum(Rmse))))


