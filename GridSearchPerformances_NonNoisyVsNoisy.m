
addpath(genpath(fullfile(pwd, 'functions')))


%% Noisy dictiony or not ? Non-noisy grid better than noisy 
% Works with both dictionaries

% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_dico = 'files/dico_grid_95760signals_4parameters.mat';
dico_size = [2000:2000:10000 20000:10000:90000];
snr_values = [inf 200 100 75 50 30 20];
repetition = 2;
signal_test_number = 10000;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(file_dico)

Rmse = zeros(length(dico_size), length(snr_values), size(Y,2), repetition); Nrmse = Rmse;

% For each iteration
for rep = 1:repetition
    
    % For each dictionary size
    for iter = 1:length(dico_size)

        fprintf('Dictionary %d / %d \t %d signals \n', iter, length(dico_size), dico_size(iter))

        rand_perm   = randperm(length(X), dico_size(iter));

        Xtrain      = abs(X(rand_perm, :));
        Ytrain      = Y(rand_perm, :);

        Xtest       = abs(X(1:signal_test_number, :));
        Ytest       = Y(1:signal_test_number, :);

        % For each snr value
%         figure
        for iter2 = 1:length(snr_values)

            Xtrain_noisy   = AddNoise(Xtrain, snr_values(iter2), 0);
            Xtest_noisy    = AddNoise(Xtest , snr_values(iter2), 0);

            Ypredict       = EstimateParametersFromGrid(Xtest_noisy, Xtrain, Ytrain);
            Ypredict_noisy = EstimateParametersFromGrid(Xtest_noisy, Xtrain_noisy, Ytrain);

            [Rmse(iter, iter2, :, rep),       Nrmse(iter, iter2, :, rep)]        = EvaluateEstimation(Ytest, Ypredict, Ytrain);
            [Rmse_noisy(iter, iter2, :, rep), Nrmse_noisy(iter, iter2, :, rep)]  = EvaluateEstimation(Ytest, Ypredict_noisy, Ytrain);
            
%             hold on
%             subplot(3, length(snr_values), iter2); plot(abs(Xtrain_noisy(1:4,:))')
%             subplot(3, length(snr_values), length(snr_values)+iter2); plot(Ytest(:,1), Ypredict(:,1), 'x')
%             subplot(3, length(snr_values), 2*length(snr_values)+iter2); plot(Ytest(:,2), Ypredict(:,2), 'x')
        end
    end
end

%% Plotting

figure
subplot(1,2,1); semilogx(dico_size, mean(mean(Nrmse, 4),3), 'o-', 'linewidth', 1.5)
hold on; semilogx(dico_size, mean(mean(Nrmse_noisy, 4),3), 'x--', 'linewidth', 1.5)
xlabel('Dictionary size'); ylabel('NRMSE'); legend(['SNR = ' + string(snr_values) 'SNR = ' + string(snr_values)]);
title('Non-noisy grid (solid lines) vs. noisy grid (dotted lines)')

subplot(1,2,2); semilogx(dico_size, mean(mean(Nrmse_noisy, 4),3) ./ mean(mean(Nrmse, 4),3), 'o-', 'linewidth', 1.5)
xlabel('Dictionary size'); ylabel('Ratio'); legend(['SNR = ' + string(snr_values) 'SNR = ' + string(snr_values)]);
title('Ratio : Noisy / Non-noisy')



