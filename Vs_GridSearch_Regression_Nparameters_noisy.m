
addpath(genpath(fullfile(pwd, 'functions')))

display = 0;


%% Regression approach Vs Grid search approach

% Load data
file_dico_grid	= 'dictionaries/multidico/2018-07-25-11:58_32-samples_5-parameters_80000-signals.mat';
file_dico_regr  = 'dictionaries/multidico/2018-07-25-09:35_32-samples_5-parameters_80000-signals.mat';

nb_signal   = 10000;
lw          = 0;
cstr.Sigma  = 'd';
cstr.Gammat = 'd';
cstr.Gammaw = '';

coord       = [1 2 3 4];

snr_values = [inf 100 20 5];

load(file_dico_regr)
Xreg        = abs(X);
Yreg        = Y(:,coord);

load(file_dico_grid)
Xgrid       = abs(X(1:length(Xreg),:));
Ygrid       = Y(1:length(Xreg),coord);
clear X Y

for rep = 1:1
    for dic = 1:1:3
        
        fprintf(['Repetition ' num2str(rep) ' - dic ' num2str(dic) '\n'])
        
        rand_perm   = randperm(length(Xreg), nb_signal);
        Xtest       = abs(Xreg(rand_perm, :));
        Ytest       = Yreg(rand_perm, :);
        
        vect    = zeros(1, length(unique(Ygrid(:,1))));
        vect(1:2^(dic-1):end) = 1;
        for c = 2:size(Ygrid,2)
            zero    = zeros(size(vect));
            vect    = repmat([vect repmat(zero, 1, 2^(dic-1) -1)], 1, floor(length(unique(Ygrid(:,c))) / (2^(dic-1)) ));
            if length(vect) ~= (length(unique(Ygrid(1,1:c), 'rows')))
                vect = [vect zeros(1, length(unique(Ygrid(:,1:c), 'rows')) - length(vect))];
            end
        end
        % If 1D % vect = repmat(vect, 1, floor(length(Ygrid) / length(unique(Ygrid(:,1)))));
        vect = repmat(vect, 1, floor(length(Ygrid) / length(vect)));
        vect = logical(vect);
        
        Xgrid_eff       = Xgrid(vect,:);
        Ygrid_eff       = Ygrid(vect,:);

        rand_perm       = randperm(length(Xreg), length(Ygrid_eff));
        Xregression_eff = Xreg(rand_perm,:);
        Yregression_eff = Yreg(rand_perm,:);        
        
        dico_len(dic)   = length(Yregression_eff);
        
        for iter = 1:length(snr_values)
            
            Xregression_eff	= AddNoise(Xregression_eff , snr_values(iter), 0);
        
            Xtest_noisy     = AddNoise(Xtest , snr_values(iter), 0);
            snr(iter,:,rep) = mean(1 ./ std(Xtest_noisy - Xtest));

            % Compute prediction
            Ypredict_grid       = EstimateParametersFromGrid(Xtest_noisy, Xgrid_eff, Ygrid_eff);
            [theta, ~]          = EstimateInverseFunction(Yregression_eff, Xregression_eff, 20, lw, 100, cstr);
            Ypredict_reg        = EstimateParametersFromModel(Xtest_noisy, theta);
            Ypredict_mean       = mean([min(Ygrid); max(Ygrid)]);

            for i = 1:size(Ygrid,2)
                Ypredict_reg(Ypredict_reg(:,i) > max(Ygrid(:,i))) = nan;
                Ypredict_reg(Ypredict_reg(:,i) < min(Ygrid(:,i))) = nan;
            end

            [Rmse_grid(dic,:,rep, iter), Nrmse_grid(dic,:,rep, iter), Mae_grid(dic,:,rep, iter)]	= EvaluateEstimation(Ytest, Ypredict_grid, Ygrid_eff);
            [Rmse_reg(dic,:,rep, iter),  Nrmse_reg(dic,:,rep, iter),  Mae_reg(dic,:,rep, iter)]     = EvaluateEstimation(Ytest, Ypredict_reg,  Yregression_eff);
            [Rmse_mean(dic,:,rep, iter), Nrmse_mean(dic,:,rep, iter), Mae_mean(dic,:,rep, iter)]	= EvaluateEstimation(Ytest, Ypredict_mean, Yregression_eff);

            if display == 1
                figure
                for i = 1:2
                    subplot(2,3,i);   plot(Ytest(:,i), Ypredict_grid(:,i), '.')
                    subplot(2,3,3+i); plot(Ytest(:,i), Ypredict_reg(:,i),  '.')
                end
            end
        end
    end
end


%%
colors = [        0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840;
    ];

figure;
subplot(221)
for i = 1:length(snr_values)
    semilogx(dico_len, squeeze(mean(Mae_grid(:,1,:,i), 3)), 'o-', 'linewidth', 1.5, 'Color', colors(i,:))
    hold on;
    semilogx(dico_len, squeeze(mean(Mae_reg(:,1,:,i), 3)), 'o--', 'linewidth', 1.5, 'Color', colors(i,:))
end
title('Coordinate 1')
ylabel('MAE'); xlabel('Dictionary sizes')


subplot(222)
for i = 1:length(snr_values)
    semilogx(dico_len, squeeze(mean(Mae_grid(:,2,:,i), 3)), 'o-', 'linewidth', 1.5, 'Color', colors(i,:))
    hold on;
    semilogx(dico_len, squeeze(mean(Mae_reg(:,2,:,i), 3)), 'o--', 'linewidth', 1.5, 'Color', colors(i,:))
end
title('Coordinate 2')
ylabel('MAE'); xlabel('Dictionary sizes')
legend({'Grid - SNR = Inf', 'Regression - SNR = Inf', 'Grid - SNR = 100', 'Regression - SNR = 100', 'Grid - SNR = 20', 'Regression - SNR = 20', 'Grid - SNR = 5', 'Regression - SNR = 5'})

subplot(223)
for i = 1:length(snr_values)
    semilogx(dico_len, squeeze(mean(Nrmse_grid(:,1,:,i), 3)), 'o-', 'linewidth', 1.5, 'Color', colors(i,:))
    hold on;
    semilogx(dico_len, squeeze(mean(Nrmse_reg(:,1,:,i), 3)), 'o--', 'linewidth', 1.5, 'Color', colors(i,:))
end
title('Coordinate 1')
ylabel('RMSE'); xlabel('Dictionary sizes')

subplot(224)
for i = 1:length(snr_values)
    semilogx(dico_len, squeeze(mean(Nrmse_grid(:,2,:,i), 3)), 'o-', 'linewidth', 1.5, 'Color', colors(i,:))
    hold on;
    semilogx(dico_len, squeeze(mean(Nrmse_reg(:,2,:,i), 3)), 'o--', 'linewidth', 1.5, 'Color', colors(i,:))
end
title('Coordinate 2')
ylabel('RMSE'); xlabel('Dictionary sizes')

%% Plotting

figure
subplot(1,2,1);
for i = 1:length(snr_values)
    semilogx(dico_len, squeeze(mean(mean(Nrmse_grid(:,:,:,i), 3),2)), 'o-', 'linewidth', 1.5, 'Color', colors(i,:))
    hold on;
    semilogx(dico_len, mean(mean(Nrmse_reg(:,:,:,i), 3),2), 'x--', 'linewidth', 1.5, 'Color', colors(i,:))
end
xlabel('Dictionary size'); ylabel('NRMSE'); 
legend({'Grid - SNR = Inf', 'Regression - SNR = Inf', 'Grid - SNR = 100', 'Regression - SNR = 100', 'Grid - SNR = 20', 'Regression - SNR = 20', 'Grid - SNR = 5', 'Regression - SNR = 5'})
title('Mean NRMSE')

subplot(1,2,2);
for i = 1:length(snr_values)
    semilogx(dico_len, squeeze(mean(mean(Nrmse_grid(:,:,:,i), 3),2) ./ mean(mean(Nrmse_reg(:,:,:,i), 3),2)), 'o-', 'linewidth', 1.5, 'Color', colors(i,:))
    hold on
end
semilogx(dico_len, [1 1 1 1], 'k', 'linewidth', 1)
xlabel('Dictionary size'); ylabel('Ratio'); legend({'SNR = Inf', 'SNR = 100','SNR = 20','SNR = 5'});
title('Ratio : Grid / Random')