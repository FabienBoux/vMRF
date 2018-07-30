% 24/01/18 
% Le but de ce script est d'étudier l'erreur d'estimation et essayer
% d'envisager une méthode afin de voir si l'on peut prévoir les erreurs
% d'estimation


%% Model learning

% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_dico       = 'files/dico_random_220000signals.mat';

signal_number   = 50000;
test_number     = 10000; 

snr_max         = 100;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(file_dico)

coord       = 1;
[~, uniq]   = unique(Y(:,coord));
Y           = Y(uniq,coord);
X           = X(uniq,:);

K           = 30;
Lt          = 0;

rand_perm   = randperm(length(X), signal_number);

Xtrain      = abs(X(rand_perm, :));
Ytrain      = Y(rand_perm, :);

rand_perm   = randperm(length(X), test_number);

Xtest       = abs(X(rand_perm, :));
Ytest       = Y(rand_perm, :);

% Create model
Ytrain_noisy = [Ytrain zeros(length(Ytrain), 1)];
Xtrain_noisy = Xtrain;
for i = 1:length(Xtrain)
    Ytrain_noisy(i,3)   = rand* snr_max;
    Xtrain_noisy(i,:)   = AddNoise(Xtrain(i,:), Ytrain_noisy(i,3), 0);
end
%%
[theta, r]	= EstimateInverseFunction(Ytrain, Xtrain, K, Lt, 30, 0);


Ypredict 	= EstimateParametersFromModel(Xtrain, theta);



%% relation between theta.Sigma and errors

err     = abs(Ytrain - Ypredict);
[~, mm] = max(r');

for i = 1:size(r,2)
    n_loc(i)   = sum(mm == i);
    err_loc(i,:) = mean(err(mm == i,:));
end

err_estim = squeeze(mean(mean(theta.Sigma,1),2));

figure
plot(err_estim, err_loc, '.')


%% Influence of noise on this error

snr     = exp(2:.5:7);
K       = 10;
Lt      = 0;

repetition = 40;

clear err err_estim err_loc n_loc
for rep = 1:repetition
    for i = 1:length(snr)
        Xtrain_noisy     = AddNoise(Xtrain, snr(i), 0);

        [theta, r]      = EstimateInverseFunction(Ytrain, Xtrain_noisy, K, Lt, 50, 0);
        Ypredict        = EstimateParametersFromModel(Xtrain_noisy, theta);

        err     = abs(Ytrain - Ypredict);
        [~, mm] = max(r');

        for j = 1:size(r,2)
            n_loc(i,j,rep)   = sum(mm == j);
            err_loc(i,j,rep) = mean(err(mm == j,:));
        end

        err_estim(i,:,rep) = squeeze(mean(mean(theta.Sigma,1),2))';
    end
end

%%
figure
subplot(121); semilogx(snr, sum(mean(err_estim,3),2), 'o');
subplot(122); semilogx(snr, sum(mean(err_loc,3),2), 'ro')
    
    
    


