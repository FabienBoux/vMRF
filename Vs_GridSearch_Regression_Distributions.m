
addpath(genpath(fullfile(pwd, 'functions')))


%% Regression approach Vs Grid search approach

% Load data
file_dico_grid        = 'files/dico_grid_216931signals.mat';
file_dico_regression  = 'files/dico_random_220000signals.mat';

nb_signal   = 1;
lw          = 0;
cstr.Sigma  = 'd';
coord       = [1 2];

snr_value   = Inf;

load(file_dico_regression)
Xreg        = abs(X);
Yreg        = Y(:,coord);

load(file_dico_grid)
Xgrid       = abs(X);
Ygrid       = Y(:,coord);
clear X Y


%% Compute regression

dic             = 5;
vect            = zeros(1, length(unique(Ygrid(:,1))));
vect(1:dic:end) = 1;

zero            = zeros(size(vect));
vect            = repmat([vect repmat(zero, 1, dic-1)], 1, ceil(length(Ygrid) / length([vect repmat(zero, 1, dic-1)])));
vect            = logical(vect);

Xgrid_eff       = Xgrid(vect,:);
Ygrid_eff       = Ygrid(vect,:);

rand_perm       = randperm(length(Xreg), length(Ygrid_eff));
Xregression_eff = Xreg(rand_perm,:);
Yregression_eff = Yreg(rand_perm,:);        

Xregression_eff = AddNoise(Xregression_eff , snr_value, 0);

[theta, ~]    	= EstimateInverseFunction(Yregression_eff, Xregression_eff, 20, lw, 100, cstr);


%% Estimation

rand_perm   = randperm(length(Xreg), nb_signal);
Xtest       = abs(Xreg(rand_perm, :));
Ytest       = Yreg(rand_perm, :);

Xtest_noisy = AddNoise(Xtest , snr_value, 0);
snr        	= 1 ./ std(Xtest_noisy - Xtest);

% Compute predictions
[Ypredict_grid, distrib_grid] = EstimateParametersFromGrid(Xtest_noisy, Xgrid_eff, Ygrid_eff);

Ypredict_reg    = EstimateParametersFromModel(Xtest_noisy, theta);

samples = 0.01:0.001:0.25;
y_samples       = [samples; 1e-6:(1e-5 - 1e-6)/(length(samples)-1):1e-5];
distrib_reg     = EstimateParametersDensityFromModel(Xtest_noisy, theta, y_samples, 0);

Ypredict_mean 	= mean([min(Ygrid); max(Ygrid)]);

for i = 1:size(Ygrid,2)
    Ypredict_reg(Ypredict_reg(:,i) > max(Ygrid(:,i))) = nan;
    Ypredict_reg(Ypredict_reg(:,i) < min(Ygrid(:,i))) = nan;
end

% Compute errors
[Rmse_grid, Nrmse_grid, Mae_grid]	= EvaluateEstimation(Ytest, Ypredict_grid, Ygrid_eff);
[Rmse_reg,  Nrmse_reg,  Mae_reg]    = EvaluateEstimation(Ytest, Ypredict_reg,  Yregression_eff);
[Rmse_mean, Nrmse_mean, Mae_mean]	= EvaluateEstimation(Ytest, Ypredict_mean, Yregression_eff);


%% Display : version 3D
figure
subplot(221);
v1 = unique(Ygrid_eff(:,1)); v2 = unique(Ygrid_eff(:,2));
surf(repmat(v1, 1, length(v2)), repmat(v2', length(v1), 1), reshape(distrib_grid, length(v1), length(v2)))
hold on;
plot3(Ypredict_grid(1), Ypredict_grid(2), 1, 'rx', 'LineWidth',3)
plot3(Ytest(1), Ytest(2), 1, 'go', 'LineWidth',3)
title(['Grid distribution : mean Nrmse = ' num2str(mean(Nrmse_grid))])
legend({'Distribution','Estimation','Real value'})

subplot(222);
[v1, v2] = meshgrid(y_samples(1,:), y_samples(2,:));
surf(v1, v2, distrib_reg)
hold on;
plot3(Ypredict_reg(1), Ypredict_reg(2), 1, 'rx', 'LineWidth',3)
plot3(Ytest(1), Ytest(2), 1, 'go', 'LineWidth',3)
title(['Regression distribution : mean Nrmse = ' num2str(mean(Nrmse_reg))])
legend({'Distribution','Estimation','Real value'})

% Display : version 2D
subplot(223);
v1 = unique(Ygrid_eff(:,1)); v2 = unique(Ygrid_eff(:,2));
imagesc(v1, v2, reshape(distrib_grid, length(v1), length(v2))')
hold on;
plot(Ypredict_grid(1), Ypredict_grid(2), 'r.', 'LineWidth', 1.5)
plot(Ytest(1), Ytest(2), 'g.', 'LineWidth', 1.5)
title(['Grid distribution : mean Nrmse = ' num2str(mean(Nrmse_grid))])
legend({'Estimation','Real value'})
set(gca,'YDir','normal')

subplot(224);
v1 = y_samples(1,:); v2 = y_samples(2,:);
imagesc(v1, v2, distrib_reg)
hold on;
plot(Ypredict_reg(1), Ypredict_reg(2), 'r.', 'LineWidth', 1.5)
plot(Ytest(1), Ytest(2), 'g.', 'LineWidth',1.5)
[~, loc1] = max(mean(distrib_reg,1)); [~, loc2] = max(mean(distrib_reg,2)); 
plot(v1(loc1), v2(loc2), 'm.', 'Linewidth', 1.5)
title(['Regression distribution : mean Nrmse = ' num2str(mean(Nrmse_reg))])
legend({'Estimation','Real value'})
set(gca,'YDir','normal')

