
addpath(genpath(fullfile(pwd, 'functions')))


%% Noisy dictiony or not ? 

% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_dico = 'files/dico_random_110224signals.mat';
dico_size = [1000:2000:10000 20000:10000:90000];
snr_values = [inf 100 50 20 10];
repetition = 1;
signal_test_number = 10000;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(file_dico)

Rmse = zeros(length(dico_size), length(snr_values), 2, repetition); Nrmse = Rmse;

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
        for iter2 = 1:length(snr_values)

            Xtrain_noisy  = AddNoise(Xtrain, 100, 0);
            Xtest         = AddNoise(Xtest , snr_values(iter2), 0);
            
            theta         = EstimateInverseFunction(Ytrain, Xtrain, 20, 0, 100, 0);
            Ypredict      = EstimateParametersFromModel(Xtest, theta);
            theta         = EstimateInverseFunction(Ytrain, Xtrain_noisy, 20, 0, 100, 0);
            Ypredict_noisy = EstimateParametersFromModel(Xtest, theta);

            [Rmse(iter, iter2, :, rep),       Nrmse(iter, iter2, :, rep)]        = EvaluateEstimation(Ytest, Ypredict, Ytrain);
            [Rmse_noisy(iter, iter2, :, rep), Nrmse_noisy(iter, iter2, :, rep)]  = EvaluateEstimation(Ytest, Ypredict_noisy(:,1:2), Ytrain);
            
            if snr_values(iter2) == 200
                figure;
                subplot(2,3,1); plot(Ytrain(:,1), Ytrain(:,2), 'x')
                ylabel('Grid')
                subplot(2,3,2); plot(Ytest(:,1), Ypredict(:,1), 'x');
                subplot(2,3,3); plot(Ytest(:,2), Ypredict(:,2), 'x')

                subplot(2,3,4); plot(Ytrain(:,1), Ytrain(:,2), 'x')
                ylabel('Random')
                subplot(2,3,5); plot(Ytest(:,1), Ypredict_noisy(:,1), 'x');
                subplot(2,3,6); plot(Ytest(:,2), Ypredict_noisy(:,2), 'x')
            end

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



