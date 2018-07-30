
addpath(genpath(fullfile(pwd, 'functions')))

%% Regression approach Vs Grid search approach

% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_dico_grid          = 'files/04-Apr-2018_64-samples_4-parameters_185367-signals.mat';
file_dico_regression    = 'files/06-Apr-2018_64-samples_4-parameters_185367-signals.mat';
% file_dico_grid          = 'files/dico_grid_216931signals.mat';
% file_dico_regression    = 'files/dico_random_220000signals.mat';

snr_values = [inf 100 20];
signal_test_number = 10000;
sampling = 5;

repetition = 1;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(file_dico_regression)
Xregression = abs(X);
Yregression = Y(:,1:2);

load(file_dico_grid)
Xgrid = abs(X);
Ygrid = Y(:,1:2);
clear X Y

Rmse_grid = zeros(sampling, length(snr_values), 2, repetition); Rmse_regression = Rmse_grid;
Nrmse_grid = Rmse_grid; Nrmse_regression = Rmse_regression;

% For each iteration
for rep = 1:repetition
    
    % For each dictionary size
    for iter = 1:sampling

        fprintf('%d / %d - Dictionary %d / %d ', rep, repetition, iter, sampling)
        
        vect    = zeros(1, length(unique(Ygrid(:,1)))); vect(1:2^(iter-1):end) = 1; %vect(1:iter:end) = 1;
        zero    = zeros(size(vect));
        vect    = repmat([vect repmat(zero, 1, 2^(iter-1) -1)], 1, ceil(length(unique(Ygrid(:,2)))));
        vect = repmat(vect, 1, ceil(length(Xgrid) / length(vect)));
        
        Xtrain_grid         = abs(Xgrid(vect(1:size(Xgrid,1)) == 1, :));
        Ytrain_grid         = Ygrid(vect(1:size(Xgrid,1)) == 1, :);

        dico_size(iter)     = size(Ytrain_grid,1);
        fprintf('\t - %d signals\n', dico_size(iter))
        
        rand_perm           = randperm(length(Xregression), dico_size(iter));
        Xtrain_reg          = abs(Xregression(rand_perm, :));
        Ytrain_reg          = Yregression(rand_perm, :);
        
        rand_perm           = randperm(length(Xregression), signal_test_number);
        Xtest               = abs(Xregression(rand_perm, :));
        Ytest               = Yregression(rand_perm, :);
        
        % For each snr value
        for iter2 = 1:length(snr_values)
            
            fprintf('\t - Snr = %d\n', snr_values(iter2))
            
            Xtest_noisy         = AddNoise(Xtest , snr_values(iter2), 0);
            
            cstr.Sigma = 'd';
            Ypredict_grid       = EstimateParametersFromGrid(Xtest_noisy, Xtrain_grid, Ytrain_grid);
            [theta, ~]          = EstimateInverseFunction(Ytrain_reg, Xtrain_reg, 20, 0, 100, cstr, 0);
            Ypredict_regression = EstimateParametersFromModel(Xtest_noisy, theta);
            
            for i = 1:size(Ygrid,2)
                Ypredict_regression(Ypredict_regression(:,i) > max(Ygrid(:,i))) = nan;
                Ypredict_regression(Ypredict_regression(:,i) < min(Ygrid(:,i))) = nan;
            end
            
            [Rmse_grid(iter, iter2, :, rep),       Nrmse_grid(iter, iter2, :, rep)]       = EvaluateEstimation(Ytest, Ypredict_grid,       Ytrain_grid);
            [Rmse_regression(iter, iter2, :, rep), Nrmse_regression(iter, iter2, :, rep)] = EvaluateEstimation(Ytest, Ypredict_regression, Ytrain_reg);
        end
    end
end


%% Plotting

colors = [           0    0.4470    0.7410
                0.8500    0.3250    0.0980
                0.9290    0.6940    0.1250
                0.4940    0.1840    0.5560
                0.4660    0.6740    0.1880
                0.3010    0.7450    0.9330
                0.6350    0.0780    0.1840];
set(groot,'defaultAxesColorOrder', colors(1:size(Rmse_grid,2),:)); clear colors

figure
subplot(2,2,1); semilogx(dico_size, mean(Rmse_grid(:,:,1), 4), 'o-', 'linewidth', 1.5)
hold on; semilogx(dico_size, mean(Rmse_regression(:,:,1), 4), 'x--', 'linewidth', 1.5)
xlabel('Dictionary size'); ylabel('RMSE (Coordinate 1)'); legend(['SNR = ' + string(snr_values) 'SNR = ' + string(snr_values)]);
title('Grid search approach (solid lines) vs. Regresion approach (dotted lines)')

subplot(2,2,2); semilogx(dico_size, mean(Rmse_grid(:,:,2), 4), 'o-', 'linewidth', 1.5)
hold on; semilogx(dico_size, mean(Rmse_regression(:,:,2), 4), 'x--', 'linewidth', 1.5)
xlabel('Dictionary size'); ylabel('RMSE (Coordinate 2)'); legend(['SNR = ' + string(snr_values) 'SNR = ' + string(snr_values)]);
title('Regular grid (solid lines) vs. Random dictionary (dotted lines)')

subplot(2,2,3); semilogx(dico_size, mean(mean(Nrmse_grid, 4),3), 'o-', 'linewidth', 1.5)
hold on; semilogx(dico_size, mean(mean(Nrmse_regression, 4),3), 'x--', 'linewidth', 1.5)
xlabel('Dictionary size'); ylabel('NRMSE'); legend(['SNR = ' + string(snr_values) 'SNR = ' + string(snr_values)]);
title('Grid search approach (solid lines) vs. Regresion approach (dotted lines)')

subplot(2,2,4); semilogx(dico_size, mean(mean(Nrmse_grid, 4),3) ./ mean(mean(Nrmse_regression, 4),3), 'o-', 'linewidth', 1.5)
xlabel('Dictionary size'); ylabel('Ratio'); legend(['SNR = ' + string(snr_values) 'SNR = ' + string(snr_values)]);
title('Ratio : Regular / Random')