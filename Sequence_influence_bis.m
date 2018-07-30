
addpath(genpath(fullfile(pwd, 'functions')))

% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_grid   	= 'files/multidico/2018-07-25-13:41_complex_dico.mat';
file_regression	= 'files/multidico/2018-07-25-11:22_complex_dico.mat';

cstr.Sigma      = '';
Lw              = 1;

coord           = [1 2 3];
seq_sizes       = [32 32 30];
snr_values      = 200;
signal_test     = 100;
nb_tests        = 100;
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
    
    
    % Select grid signals
    Xgrid               = Xdico(:,nmat(seq,:));
    Ygrid               = Ydico;
    
    % Select trainning data and train model
    Xtrain              = Xregr(:,nmat(seq,:));
    Ytrain              = Yregr;
    [theta, ~]          = EstimateInverseFunction(Ytrain, Xtrain, 20, Lw, 100, cstr, 0);
    
    for rep = 1:nb_tests
        % Generate test signals (adding noise)
        rand_perm           = randperm(length(Xregr), signal_test);
        Xtest               = Xregr(rand_perm, nmat(seq,:));
        Ytest               = Yregr(rand_perm, :);
        Xtest               = AddNoise(Xtest , snr_values, 0);

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

%%




