
addpath(genpath(fullfile(pwd, 'functions')))

%% Random or regular grid ? 

% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_dico_regular = 'files/dico_grid_216931signals.mat';
file_dico_random  = 'files/dico_random_110224signals.mat';

snr_values = [inf 200 100 75 50 30 20];
signal_test_number = 10000;
sampling= 6;

repetition = 1;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(file_dico_random)
Xrandom = X;
Yrandom = Y;
load(file_dico_regular)

Rmse = zeros(sampling, length(snr_values), 2, repetition); Nrmse = Rmse;

% For each iteration
for rep = 1:repetition
    
    % For each dictionary size
    for iter = 1:sampling

        fprintf('Dictionary %d / %d \n', iter, sampling)
        
        vect = zeros(1, length(unique(Y(:,1)))); vect(1:2^iter:end) = 1; %vect(1:iter:end) = 1;
        zero = zeros(size(vect));
        vect = repmat([vect repmat(zero, 1, 2^iter -1)], 1, ceil(length(unique(Y(:,2)))));
        
        Xtrain = abs(X(vect(1:size(X,1)) == 1, :));
        Ytrain = Y(vect(1:size(X,1)) == 1, :);

        dico_size(iter) = size(Ytrain,1);
        rand_perm       = randperm(length(Xrandom), dico_size(iter));
        Xtrain_random   = abs(Xrandom(rand_perm, :));
        Ytrain_random   = Yrandom(rand_perm, :);
        
        rand_perm       = randperm(length(Xrandom), signal_test_number);
        Xtest           = abs(Xrandom(rand_perm, :));
        Ytest           = Yrandom(rand_perm, :);
        
        % For each snr value
        for iter2 = 1:length(snr_values)
            Xtest_noisy     = AddNoise(Xtest , snr_values(iter2), 0);

            Ypredict        = EstimateParametersFromGrid(Xtest_noisy, Xtrain, Ytrain);
            Ypredict_random = EstimateParametersFromGrid(Xtest_noisy, Xtrain_random, Ytrain_random);

            [Rmse(iter, iter2, :, rep),         Nrmse(iter, iter2, :, rep)]         = EvaluateEstimation(Ytest, Ypredict, Ytrain);
            [Rmse_random(iter, iter2, :, rep),  Nrmse_random(iter, iter2, :, rep)]  = EvaluateEstimation(Ytest, Ypredict_random, Ytrain);
            
            
%             if snr_values(iter2) == 100
%                 figure;
%                 subplot(2,3,1); plot(Ytrain(:,1), Ytrain(:,2), 'x')
%                 hold on
%                 plot(Ypredict(1:10,1), Ypredict(1:10,2), 'go')
%                 plot(Ytest(1:10,1), Ytest(1:10,2), 'ro')
%                 ylabel('Grid')
%                 subplot(2,3,2); plot(Ytest(:,1), Ypredict(:,1), 'x'); xlim([0 0.25]); ylim([0 0.25])
%                 subplot(2,3,3); plot(Ytest(:,2), Ypredict(:,2), 'x'); xlim([1 10]); ylim([1 10])
% 
%                 subplot(2,3,4); plot(Ytrain_random(:,1), Ytrain_random(:,2), 'x')
%                 hold on
%                 plot(Ypredict_random(1:10,1), Ypredict_random(1:10,2), 'go')
%                 plot(Ytest(1:10,1), Ytest(1:10,2), 'ro')
%                 ylabel('Random')
%                 subplot(2,3,5); plot(Ytest(:,1), Ypredict_random(:,1), 'x'); xlim([0 0.25]); ylim([0 0.25])
%                 subplot(2,3,6); plot(Ytest(:,2), Ypredict_random(:,2), 'x'); xlim([1 10]); ylim([1 10])
%             end
        end
    end
end

%% Plotting

figure
subplot(2,2,1); semilogx(dico_size, mean(Rmse(:,:,1), 4), 'o-', 'linewidth', 1.5)
hold on; semilogx(dico_size, mean(Rmse_random(:,:,1), 4), 'x--', 'linewidth', 1.5)
xlabel('Dictionary size'); ylabel('RMSE (Coordinate 1)'); legend(['SNR = ' + string(snr_values) 'SNR = ' + string(snr_values)]);
title('Regular grid (solid lines) vs. Random dictionary (dotted lines)')
subplot(2,2,2); semilogx(dico_size, mean(Rmse(:,:,2), 4), 'o-', 'linewidth', 1.5)
hold on; semilogx(dico_size, mean(Rmse_random(:,:,2), 4), 'x--', 'linewidth', 1.5)
xlabel('Dictionary size'); ylabel('RMSE (Coordinate 2)'); legend(['SNR = ' + string(snr_values) 'SNR = ' + string(snr_values)]);
title('Regular grid (solid lines) vs. Random dictionary (dotted lines)')

subplot(2,2,3); semilogx(dico_size, mean(mean(Nrmse, 4),3), 'o-', 'linewidth', 1.5)
hold on; semilogx(dico_size, mean(mean(Nrmse_random, 4),3), 'x--', 'linewidth', 1.5)
xlabel('Dictionary size'); ylabel('NRMSE'); legend(['SNR = ' + string(snr_values) 'SNR = ' + string(snr_values)]);
title('Regular grid (solid lines) vs. Random dictionary (dotted lines)')
subplot(2,2,4); semilogx(dico_size, mean(mean(Nrmse, 4),3) ./ mean(mean(Nrmse_random, 4),3), 'o-', 'linewidth', 1.5)
xlabel('Dictionary size'); ylabel('Ratio'); legend(['SNR = ' + string(snr_values) 'SNR = ' + string(snr_values)]);
title('Ratio : Regular / Random')
