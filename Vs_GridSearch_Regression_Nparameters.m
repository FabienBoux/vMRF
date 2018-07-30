
addpath(genpath(fullfile(pwd, 'functions')))

display = 0;


%% Regression approach Vs Grid search approach

% Load data
file_dico_grid        = 'files/17-07-2018_dico_grid_4params_160000signals.mat';
file_dico_regression  = 'files/16-07-2018_dico_rand_4params_160000signals.mat';

nb_signal   = 10000;
lw          = 1;
cstr.Sigma  = 'd';
noise       = 100;

load(file_dico_regression)
Xregr   = abs(X);
Yregr   = Y;

load(file_dico_grid)
Xgrid  	= abs(X);
Ygrid  	= Y;
clear X Y

for rep = 1:3
    for dic = 1:1:3
        
        fprintf(['Repetition ' num2str(rep) ' - dic ' num2str(dic) '\n'])
        
        rand_perm   = randperm(length(Xregr), nb_signal);
        Xtest       = abs(Xregr(rand_perm, :));
        Ytest       = Yregr(rand_perm, :);
        
        Xtest_noisy = AddNoise(Xtest , noise, 0);
        snr(dic,:,rep) = 1 ./ std(Xtest_noisy - Xtest);
        
        vect = [];
        for coord = 1:size(Ygrid,2)
            s       = length(unique(Ygrid(:,coord)));
            if isempty(vect)
                vect = zeros(1,s);
                vect(1:dic:end) = 1;
            else
                tmp     = [vect repmat(zeros(1,length(vect)), 1,dic-1)];
                tmp     = repmat(tmp, 1,ceil(s/dic));
                vect    = zeros(1,length(vect)*s);
                vect	= tmp(1:length(vect));
            end
        end
        vect = logical(vect);
        
        Xgrid_eff       = Xgrid(vect,:);
        Ygrid_eff       = Ygrid(vect,:);
        
        if length(Ygrid_eff) >= length(Xregr)
            rand_perm       = randperm(length(Xregr), length(Xregr));
        else
            rand_perm       = randperm(length(Xregr), length(Ygrid_eff));
        end
        
        Xregr_eff       = AddNoise(Xregr(rand_perm,:), noise, 0);
        Yregr_eff       = Yregr(rand_perm,:);        
        
        dico_len(dic)   = length(Yregr_eff);

        % Compute prediction
        Ypredict_grid       = EstimateParametersFromGrid(Xtest_noisy, Xgrid_eff, Ygrid_eff);
        
        [theta, ~]          = EstimateInverseFunction(Yregr_eff, Xregr_eff, 20, lw, 200, cstr);
        Ypredict_regr       = EstimateParametersFromModel(Xtest_noisy, theta);
        Ypredict_rand       = min(Ygrid) + rand(length(Ypredict_regr),size(Ygrid,2)) .* (max(Ygrid) - min(Ygrid));
        Ypredict_mean       = mean([min(Ygrid); max(Ygrid)]);
        
        for i = 1:size(Ygrid,2)
            Ypredict_regr(Ypredict_regr(:,i) > max(Ygrid(:,i))) = nan;
            Ypredict_regr(Ypredict_regr(:,i) < min(Ygrid(:,i))) = nan;
        end
        
        [Rmse_grid(dic,:,rep), Nrmse_grid(dic,:,rep)]	= EvaluateEstimation(Ytest, Ypredict_grid, Ygrid_eff);
        [Rmse_regr(dic,:,rep), Nrmse_regr(dic,:,rep)]   = EvaluateEstimation(Ytest, Ypredict_regr,  Yregr_eff);
        [Rmse_rand(dic,:,rep), Nrmse_rand(dic,:,rep)]   = EvaluateEstimation(Ytest, Ypredict_rand, Yregr_eff);
        [Rmse_mean(dic,:,rep), Nrmse_mean(dic,:,rep)]	= EvaluateEstimation(Ytest, Ypredict_mean, Yregr_eff);
        
        if display == 1
            figure
            for i = 1:2
                subplot(2,3,i);   plot(Ytest(:,i), Ypredict_grid(:,i), '.')
                subplot(2,3,3+i); plot(Ytest(:,i), Ypredict_regr(:,i), '.')
            end
        end
    end
end


%%
figure
for coord = 1:size(Rmse_grid,2)
    subplot(2,size(Rmse_grid,2),coord); semilogx(dico_len, mean(Rmse_grid(:,coord,:),3), 'o-')
    hold on; semilogx(dico_len, mean(Rmse_regr(:,coord,:),3), 'o-')
    semilogx(dico_len, mean(Rmse_rand(:,coord,:),3), 'o-')
    semilogx(dico_len, mean(Rmse_mean(:,coord,:),3), 'o-')
    legend({'Grid', 'Regression', 'RandomUnif', 'Mean'})
    title(['Coord ' num2str(coord)])
    if coord == 1, ylabel('Rmse'); end
    subplot(2,size(Rmse_grid,2),size(Rmse_grid,2)+coord); semilogx(dico_len, mean(Nrmse_grid(:,coord,:),3), 'o-')
    hold on; semilogx(dico_len, mean(Nrmse_regr(:,coord,:),3), 'o-')
    semilogx(dico_len, mean(Nrmse_rand(:,coord,:),3), 'o-')
    semilogx(dico_len, mean(Nrmse_mean(:,coord,:),3), 'o-')
    legend({'Grid', 'Regression', 'RandomUnif', 'Mean'})
    ylim([0 1.5])
    if coord == 1, ylabel('Nrmse'); end
end


%%
figure
semilogx(dico_len, mean(mean(Nrmse_grid,3),2), 'o-', 'LineWidth', 1.5)
hold on
semilogx(dico_len, mean(mean(Nrmse_regr,3),2), 'o-', 'LineWidth', 1.5)
semilogx(dico_len, mean(mean(Nrmse_rand,3),2), 'o-', 'LineWidth', 1.5)
semilogx(dico_len, mean(mean(Nrmse_mean,3),2), 'o-', 'LineWidth', 1.5)
legend({'Grid', 'Regression', 'RandomUnif', 'Mean'})
% ylim([0 1.5])


