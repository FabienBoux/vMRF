% Note lines ~30/31 used because of the dictionary, else remove

%%

file_grid = 'dictionaries/multidico/18-08-21/2018-08-22-10:08_32-samples_4-parameters_36000-signals.mat';
file_rand = 'dictionaries/multidico/18-08-21/2018-08-22-09:53_32-samples_4-parameters_36000-signals.mat';

subsamplmax = 6;

K           = 20;
Lw          = 1;
cstr.Sigma  = 'd';
cstr.Gammat = 'd';
cstr.Gammaw = '';

snr         = [20 200];
snr_train   = snr;

%%


load(file_rand)
Xtest   = X;
Xtest   = AddNoise(Xtest, snr(1)+(snr(2)-snr(1))*rand(size(Xtest,1),1), 0);
Ytest   = Y;

Xtrain  = Xtest;
Ytrain  = Ytest;

% Xtest   = Xtest(Ytest(:,1) <=0.2,:);
% Ytest   = Ytest(Ytest(:,1) <=0.2,:);

load(file_grid)


% figure
for j = 1:subsamplmax  %j = 1 corresponding to all dico, j = 2 the half, ...
    
    Xgrid = X;
    Ygrid = Y;
    
    % Reduce grid depending on j
    vect = true(size(Ygrid,1),1);
    for c = 1:size(Ygrid,2)
        tmp_vect{j,c} = false(size(Ygrid,1),1);
        val     = unique(Ygrid(:,c));
        val     = val(1:j:end);
        if length(val) < 3
            tmp_vect{j,c} = tmp_vect{j-1,c};
        else
            for v = 1:length(val), tmp_vect{j,c} = tmp_vect{j,c} | (Ygrid(:,c) == val(v)); end
        end
        vect = vect & tmp_vect{j,c};
    end
    Xgrid = Xgrid(vect,:);
    Ygrid = Ygrid(vect,:);
    
    dico_size(j) = size(Xgrid,1);
    for c = 1:size(Ygrid,2), step_values(j,c) = mean(unique((diff(unique(Ygrid(:,c)))))); end
    
    % Compute grid search and prediciton
    Ypredict_grid   = EstimateParametersFromGrid(Xtest, Xgrid, Ygrid);
    
    
    %Xtrain	= AddNoise(Xtrain, snr_train(1)+(snr_train(2)-snr_train(1))*rand(size(Xtrain,1),1), 0);
    
    % Reduce learning dataset
    vect    = randperm(size(Xtrain,1), dico_size(j));
    Xtrain  = Xtrain(vect,:);
    Ytrain  = Ytrain(vect,:);
    
    Xtrain  = Xgrid;
    Ytrain  = Ygrid;
    
    % Compute regression and prediction
    [theta,~] 	= EstimateInverseFunction(Ytrain, Xtrain, K, Lw, 200, cstr, 0);
    Ypredict_regr	= EstimateParametersFromModel(Xtest, theta, 0);

    for i = 1:size(Ytrain,2)
        Ypredict_regr(Ypredict_regr(:,i) > max(Ytrain(:,i)),i) = nan;
        Ypredict_regr(Ypredict_regr(:,i) < min(Ytrain(:,i)),i) = nan;
    end
    
    % Evaluate predictions
    [Rmse_grid(j,:), Nrmse_grid(j,:), Mae_grid(j,:)] = EvaluateEstimation(Ytest, Ypredict_grid, Ygrid);
    [Rmse_regr(j,:), Nrmse_regr(j,:), Mae_regr(j,:)] = EvaluateEstimation(Ytest, Ypredict_regr, Ytrain);
    
    % Temporary display : use it to illustrate the subsampling of the grid
%     subplot(4,subsamplmax,j); plot(Ygrid(:,1), Ygrid(:,2), 'o')
%     subplot(4,subsamplmax,subsamplmax+j); plot(Ygrid(:,2), Ygrid(:,3), 'o')
%     
%     subplot(4,subsamplmax,2*subsamplmax+j); plot(Ygrid(:,1), Ygrid(:,2), 'o')
%     subplot(4,subsamplmax,3*subsamplmax+j); plot(Ygrid(:,2), Ygrid(:,3), 'o')
end


%% 

figure
for c = 1:size(Y,2)
    subplot(2,2,c)
    plot(step_values(:,c), Mae_grid(:,c), 'o-', 'LineWidth',1.5)
    hold on 
    mdl = fitlm(step_values(:,c), Mae_grid(:,c));
    p = mdl.Coefficients.Estimate;
    plot(step_values(:,c), step_values(:,c)*p(2) + p(1), 'LineWidth',1.5)
    title(['error = ' num2str(p(2),2) ' step + ' num2str(p(1),2) ' / R² = ' num2str(mdl.Rsquared.Ordinary,3)])
end

figure
for c = 1:size(Y,2)
    subplot(2,2,c)
    plot(step_values(:,c), Mae_regr(:,c), 'o-', 'LineWidth',1.5)
    hold on 
    mdl = fitlm(step_values(:,c), Mae_regr(:,c));
    p = mdl.Coefficients.Estimate;
    plot(step_values(:,c), step_values(:,c)*p(2) + p(1), 'LineWidth',1.5)
    title(['error = ' num2str(p(2),2) ' step + ' num2str(p(1),2) ' / R² = ' num2str(mdl.Rsquared.Ordinary,3)])
end


%%
figure
for c = 1:size(Y,2)
    subplot(2,2,c)
    semilogx(dico_size, Mae_grid(:,c), 'o-')
    hold on
    semilogx(dico_size, Mae_regr(:,c), 'o-')
    legend({'Grid', 'Regression'})
    title(['MAE coord: ' num2str(c)])
end

figure
semilogx(dico_size, mean(Nrmse_grid,2), 'o-')
hold on
semilogx(dico_size, mean(Nrmse_regr,2), 'o-')
legend({'Grid', 'Regression'})
title('Average NRMSE')


