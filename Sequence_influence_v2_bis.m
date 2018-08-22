% SEQUENCE_INFLUENCE
%
% Using a dictionary composed of MGEFIDSE, MSME and MGE (pre/post) signals,
% the idea is to compare estimation (average NRMSE distributions) using  
% each sequence only or a combinaison of sequences
%
% Note: due to some errors with dictionaries (BVf interval definition), a
% line was added to remove some test signals
%
% Fabien Boux - 08/2018


addpath(genpath(fullfile(pwd, 'functions')))

% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_grid   	= 'dictionaries/multidico/18-08-22-2/2018-08-22-17:14-24_complex_dico_2.mat';
file_regression	= 'dictionaries/multidico/18-08-22-2/2018-08-22-17:05-59_complex_dico_2.mat';
seq_names       = {'GEFIDSE','MGE','MSME'};

cstr.Sigma      = 'd';
cstr.Gammat     = 'd';
cstr.Gammaw     = '';
Lw              = 0;
K               = 20;

coord           = [1 2 3];
seq_sizes       = [32 32 30];
signal_test     = 100;
nb_tests        = 1000;

snr_train       = [0 inf];
snr_test        = [0 inf];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Load dicos
load(file_regression)
Xregr = X;
Yregr = Y(:,coord);

load(file_grid)
Xdico = X;
Ydico = Y(:,coord);
clear X Y


% Convert seq_sizes to logical vectors
clear mat nmat
vect    = ones(1,length(seq_sizes));
mat     = vect;
for i = 1:length(seq_sizes)-1
    vect(1:i)   = 0;
    mat         = [mat; unique(perms(vect), 'rows')];
end
for i = 1:size(mat,1), nmat(i,:) = repelem(mat(i,:), seq_sizes); end
nmat   	= logical(repmat(nmat, 1,2));


for seq = 1:size(mat,1)
    
    fprintf(['Combination : ' num2str(seq) '/' num2str(size(mat,1)) '\n'])
    
    % Select grid signals
    Xgrid               = Xdico(:,nmat(seq,:));
    Ygrid               = Ydico;
    
    % Select trainning data and train model
    Xtrain              = Xregr(:,nmat(seq,:));
    Xtrain              = AddNoise(Xtrain, snr_train(1)+(snr_train(2)-snr_train(1))*rand(size(Xtrain,1),1), 0);
    Ytrain              = Yregr;
    [theta,~]           = EstimateInverseFunction(Ytrain, Xtrain, K, Lw, 100, cstr, 0);
    
    for rep = 1:nb_tests
        % Generate test signals (adding noise)
        rand_perm           = randperm(length(Xregr), signal_test);
        Xtest               = Xregr(rand_perm, nmat(seq,:));
        Ytest               = Yregr(rand_perm, :);
        Xtest               = AddNoise(Xtest , snr_test(1)+(snr_test(2)-snr_test(1))*rand(size(Xtest,1),1), 0); 
        
        % Prediction
        Ypredict_grid       = EstimateParametersFromGrid(Xtest, Xgrid, Ygrid);
        Ypredict_regr       = EstimateParametersFromModel(Xtest, theta, 0);

        for i = 1:size(Ygrid,2)
            Ypredict_regr(Ypredict_regr(:,i) > max(Ygrid(:,i)),i) = nan;
            Ypredict_regr(Ypredict_regr(:,i) < min(Ygrid(:,i)),i) = nan;
        end

        % Evaluate estimations
        [Rmse_grid(seq,rep,:), Nrmse_grid(seq,rep,:), Mae_grid(seq,rep,:)] = EvaluateEstimation(Ytest, Ypredict_grid, Ygrid);
        [Rmse_regr(seq,rep,:), Nrmse_regr(seq,rep,:), Mae_regr(seq,rep,:)] = EvaluateEstimation(Ytest, Ypredict_regr, Ytrain);
    end
end


%% NRMSE

figure;
for i = 1:size(Nrmse_grid,1)
    subplot(size(Nrmse_grid,1),2,2*i-1)
    hist(squeeze(Nrmse_grid(i,:,:)), 0:0.01:1.4)
    xlim([0 1.4])
    hold on; line([1 1], ylim, 'LineWidth', 1.5, 'Color', 'k');
    if i == 1
        title('NRMSE histogram - Grid search method');
        legend({'BVf','VSI','T2','Limit of interest'})
        ylabel(seq_names)
    else
        ylabel(mat(i,:))
    end
    
    
    subplot(size(Nrmse_grid,1),2,2*i)
    hist(squeeze(Nrmse_regr(i,:,:)), 0:0.01:1.4)
    xlim([0 1.4])
    hold on; line([1 1], ylim, 'LineWidth', 1.5, 'Color', 'k');
    if i == 1
        title('NRMSE histogram - Learning method');
        legend({'BVf','VSI','T2','Limit of interest'})
    end
end

%% MAE

figure;
for i = 1:size(Mae_grid,1)
    
 
    for c = 1:size(Mae_grid,3)
        subplot(size(Mae_grid,1),size(Mae_grid,3),size(Mae_grid,3)*(i-1)+c)
        
        switch c
            case 1
                b = (0:0.00015:0.015) *1e2;
                hist([squeeze(Mae_grid(i,:,c)); squeeze(Mae_regr(i,:,c))]'*1e2, b)
                xlim([b(1) b(end)])
                if i == 1, title('BVf (in %)'); end
            case 2
                b = 0:2e-8:2e-6;
                hist([squeeze(Mae_grid(i,:,c)); squeeze(Mae_regr(i,:,c))]', b)
                xlim([b(1) b(end)])
                if i == 1, title('VSI (in µm)'); end
            case 3
                b = (0:0.0002:0.03) * 1e3;
                hist([squeeze(Mae_grid(i,:,c)); squeeze(Mae_regr(i,:,c))]' *1e3, b)
                xlim([b(1) b(end)])
                if i == 1, title('T_2 (in ms)'); end
            case 4
                break
                %TODO
        end
    end
    
    if i == 1, legend({'Grid search','Regression approach'}); end
end

