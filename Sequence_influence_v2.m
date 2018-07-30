
addpath(genpath(fullfile(pwd, 'functions')))

% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coord           = [1 2 3];
ratio           = 0;

seq_names       = {'MGEFIDSE','MSME','MGE'};
                  % First = random / second = grid
seq_files       = {{'files/multidico/2018-07-25-09:35_32-samples_5-parameters_80000-signals',...
                    'files/multidico/2018-07-25-11:58_32-samples_5-parameters_80000-signals'},...
                   {'files/multidico/2018-07-25-10:05_32-samples_5-parameters_80000-signals',...
                    'files/multidico/2018-07-25-12:31_32-samples_5-parameters_80000-signals'},...
                   {'files/multidico/2018-07-25-11:22_30-samples_5-parameters_80000-signals',...
                    'files/multidico/2018-07-25-13:41_30-samples_5-parameters_80000-signals'}};
seq_sizes       = [32 32 30];

cstr.Sigma      = 'd';
Lw              = 1;

signal_test     = 100;
nb_tests        = 100;
snr_train       = [20 500];
snr_test        = [20 100];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(seq_files) ~= length(seq_names)
    error('')
end

for seq = 1:length(seq_files)
    
    fprintf(['Sequence: ' seq_names{seq} '\n'])
    
    % Load grid signals
    load(seq_files{seq}{2})
    if ratio == 1, X = X(:,1:end/2) ./ X(:,end/2+1:end); end
    Xgrid               = X;
    Ygrid               = Y(:,coord);
    clear X Y
    
    % Load training set + train model
    load(seq_files{seq}{1})
    if ratio == 1, X = X(:,1:end/2) ./ X(:,end/2+1:end); end
    Xtrain              = X;
    Xtrain              = [Xtrain; AddNoise(Xtrain, snr_train(1)+(snr_train(2)-snr_train(1))*rand(size(Xtrain,1),1), 0)];
    Ytrain              = [Y(:,coord); Y(:,coord)];
    [theta, ~]          = EstimateInverseFunction(Ytrain, Xtrain, 20, Lw, 100, cstr, 0);
   
    for rep = 1:nb_tests
        % Generate test signals
        rand_perm           = randperm(length(X), signal_test);
        Xtest               = X(rand_perm,:);
        Ytest               = Y(rand_perm,coord);
        Xtest               = AddNoise(Xtest, snr_test(1)+(snr_test(2)-snr_train(1))*rand(size(Xtest,1),1), 0);

        % Predict using regression
        Ypredict_regr       = EstimateParametersFromModel(Xtest, theta, 0);
        for i = 1:size(Y(:,coord),2)
            Ypredict_regr(Ypredict_regr(:,i) > max(Ygrid(:,i)),i) = nan;
            Ypredict_regr(Ypredict_regr(:,i) < min(Ygrid(:,i)),i) = nan;
        end

        % Predict using grid
        Ypredict_grid       = EstimateParametersFromGrid(Xtest, Xgrid, Ygrid);
        
        % Evaluate estimations
        [Rmse_grid(seq,rep,:), Nrmse_grid(seq,rep,:), Mae_grid(seq,rep,:)] = EvaluateEstimation(Ytest, Ypredict_grid, Ygrid);
        [Rmse_regr(seq,rep,:), Nrmse_regr(seq,rep,:), Mae_regr(seq,rep,:)] = EvaluateEstimation(Ytest, Ypredict_regr, Ytrain);
    end
    clear X Y
end

%%

for i = 1:size(Nrmse_grid,1)
    fprintf('NRMSE of sequence %s :\n', seq_names{i})
    fprintf(['\t using grid search method:\t' repmat('%.2f ',1,size(Nrmse_grid,3)) ' +/- ' repmat('%.3f ',1,size(Nrmse_grid,3)) '\n'], mean(squeeze(Nrmse_grid(i,:,:))), std(squeeze(Nrmse_grid(i,:,:))));
    fprintf(['\t using regression method:\t'  repmat('%.2f ',1,size(Nrmse_regr,3)) ' +/- ' repmat('%.3f ',1,size(Nrmse_regr,3)) '\n'], mean(squeeze(Nrmse_regr(i,:,:))), std(squeeze(Nrmse_regr(i,:,:))));
end

