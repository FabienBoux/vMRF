
addpath(genpath(fullfile(pwd, 'functions')))

%% Regression approach Vs Grid search approach

% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dico_regression     = 'dictionaries/multidico/2018-07-23-09:47_rand_91-samples_4-parameters_20000-signals.mat';
dico_grid           = 'dictionaries/multidico/2018-07-23-09:22_grid_91-samples_4-parameters_20000-signals.mat';

repetition          = 1;
signal_test_number  = 200;
snr_lvl             = [50 200];
cstr.Sigma          = 'd';
Lw                  = 0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(dico_regression)
Xregr   = X;
Yregr   = Y;

load(dico_grid)
Xgrid   = X;
Ygrid   = Y;
clear   X Y


% For each parameter size
for rep = 1:repetition
    
    fprintf([num2str(rep) 'th repetition\n'])
    
    parfor_progress(size(Yregr,2));
    parfor iter = 1:size(Yregr,2)

        % Build test dataset
        rand_perm           = randperm(length(Xregr), signal_test_number);
        Xtest{iter}         = abs(Xregr(rand_perm, :));
        Ytest{iter}         = Yregr(rand_perm, 1:iter);
        
        %[Xtest{iter}, snr_test{iter}] = AddNoise(Xtest{iter}, snr_lvl(1) + ones(size(Xtest{iter},1),1)*(snr_lvl(2) - snr_lvl(1)), 0);
        
        % Build learning dataset
        vect                = ones(size(Xregr,1),1); vect(rand_perm) = 0;
        Xtrain{iter}        = Xregr(logical(vect), :);
%         for i = 1:size(Xtrain{iter},1)
%         	[Xtrain{iter}, snr_train{iter}] = AddNoise(Xtrain{iter}, snr(1) + ones(size(Xtrain{iter},1),1)*(snr_lvl(2) - snr_lvl(1)), 0);
%         end
        Ytrain{iter}        = Yregr(logical(vect), 1:iter);
        dico_size(iter)     = size(Ytrain{iter},1);

        % Train model
        [theta{iter}, ~]  	= EstimateInverseFunction(Ytrain{iter}, Xtrain{iter}, 20, Lw, 100, cstr, 0);

        % Predict using the trained model
        Ypredict_regr{iter} = EstimateParametersFromModel(Xtest{iter}, theta{iter}); 
        for i = 1:iter
            Ypredict_regr{iter}(Ypredict_regr{iter}(:,i) > max(Ygrid(:,i))) = nan;
            Ypredict_regr{iter}(Ypredict_regr{iter}(:,i) < min(Ygrid(:,i))) = nan;
        end
%         for i = iter+1:5
%             Ypredict_regr{iter} = nan(size(Ypredict_regr{iter}));
%         end

        % Predict using the grid
        Ypredict_grid{iter}      = EstimateParametersFromGrid(Xtest{iter}, Xgrid, Ygrid(:,1:iter), 0);

        % Compute estimation error
        [Rmse_regr{iter,rep}, Nrmse_regr{iter,rep}, Mae_regr{iter,rep}] = EvaluateEstimation(Ytest{iter}, Ypredict_regr{iter}, Ytrain{iter});
        [Rmse_grid{iter,rep}, Nrmse_grid{iter,rep}, Mae_grid{iter,rep}] = EvaluateEstimation(Ytest{iter}, Ypredict_grid{iter}, Ytrain{iter});

        parfor_progress;
    end
end
parfor_progress(0);



%% disp

figure
clear temp
for iter = 1:size(Ygrid,2)
    tmp = [];
    tmp_regr = [];
    for rep = 1:repetition
        tmp = [tmp; Nrmse_grid{iter,rep}];
        tmp_regr = [tmp_regr; Nrmse_regr{iter,rep}];
    end
    temp{iter} = zeros(size([tmp tmp_regr]));
    temp{iter}(:,1:2:end) = tmp;
    temp{iter}(:,2:2:end) = tmp_regr;
    
    subplot(ceil(size(Ygrid,2)/2),ceil(size(Ygrid,2)/2),iter);
    boxplot(temp{iter}, 'Labels', repmat({'Grid search','Regression'},1,iter))
    title(['NRMSE for simultaneous estimation of ' num2str(iter) ' parameters'])
end

%%
figure
clear temp
for iter = 1:size(Ygrid,2)
    tmp = [];
    for rep = 1:repetition
        tmp = [tmp; mean(Nrmse_grid{iter,rep}) - mean(Nrmse_regr{iter,rep})];
    end
    temp(:,iter) = tmp;
end
errorbar(mean(temp), std(temp), 'o-', 'LineWidth', 1.5)
title(['Mean of NRMSE_{Grid search} - NRMSE_{Regression approach} on ' num2str(repetition) ' repetitions'])
xlabel('Parameters')


