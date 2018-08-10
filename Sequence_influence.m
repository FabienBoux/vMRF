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
file_grid   	= 'dictionaries/multidico/2018-07-25-13:41_complex_dico.mat';
file_regression	= 'dictionaries/multidico/2018-07-25-11:22_complex_dico.mat';

cstr.Sigma      = 'd';
Lw              = 0;

coord           = [1 2 3 4];
seq_sizes       = [32 32 30];
snr_values      = Inf;
signal_test     = 1000;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load dico
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
    
    % Generate test signals (adding noise)
    rand_perm           = randperm(length(Xregr), signal_test);
    Xtest               = Xregr(rand_perm, nmat(seq,:));
    Ytest               = Yregr(rand_perm, :);
    %%%%%%%% lines addded due to the specific problem %%%%%%%%%%%%%%%%%%%%%
    remove_samples      = Ytest(:,1) <= 0.2;
    Ytest               = Ytest(remove_samples,:);
    Xtest               = Xtest(remove_samples,:);    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Xtest               = AddNoise(Xtest , snr_values, 0);
    
    % Predict using grid
    Xgrid               = Xdico(:,nmat(seq,:));
    Ygrid               = Ydico;
    Ypredict_grid       = EstimateParametersFromGrid(Xtest, Xgrid, Ygrid);
    
    % Predict using regression
    Xtrain              = Xregr(:,nmat(seq,:));
    Ytrain              = Yregr;
    [theta, ~]          = EstimateInverseFunction(Ytrain, Xtrain, 20, Lw, 100, cstr, 0);
    Ypredict_regr       = EstimateParametersFromModel(Xtest, theta, 0);
    
    for i = 1:size(Ygrid,2)
        Ypredict_regr(Ypredict_regr(:,i) > max(Ygrid(:,i)),i) = nan;
        Ypredict_regr(Ypredict_regr(:,i) < min(Ygrid(:,i)),i) = nan;
    end
    
    % Evaluate estimations
    [Rmse_grid(seq,:), Nrmse_grid(seq,:), Mae_grid(seq,:)] = EvaluateEstimation(Ytest, Ypredict_grid, Ygrid);
    [Rmse_regr(seq,:), Nrmse_regr(seq,:), Mae_regr(seq,:)] = EvaluateEstimation(Ytest, Ypredict_regr, Ytrain);
end





