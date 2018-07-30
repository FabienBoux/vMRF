% 22/01/18 
% Le but est de savoir si un paramètre de simulation peut être
% automatiquement appris de manière non supervisé (en précisant une 
% variable latente)


%% First try with BVf and VSI

% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_dico       = 'files/dico_random_220000signals.mat';

signal_number   = 50000;
test_number     = 20000; 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(file_dico)

K           = 15;
Lt          = 1;

rand_perm   = randperm(length(X), signal_number);

Xtrain      = abs(X(rand_perm, :));
Ytrain      = Y(rand_perm, :);

rand_perm   = randperm(length(X), test_number);

Xtest       = abs(X(rand_perm, :));
Ytest       = Y(rand_perm, :);

[theta, ~]  = EstimateInverseFunction(Ytrain(:,2), Xtrain, K, Lt, 30, 0);

Ypredict    = EstimateParametersFromModel(Xtest, theta);


%% Plot
Ytest = Ytest(:,[2 1]);

C = corr(Ypredict, Ytest);

figure
subplot(121); plot(Ytest(:,1), Ypredict(:,1), 'x')
title('Corr = ' +  string(C(1,1)));
xlabel('Real values'); ylabel('Prediction')

subplot(122); plot(Ytest(:,2), Ypredict(:,2), 'x')
title('Corr = ' +  string(C(2,2)));
xlabel('Real values'); ylabel('Prediction')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second try with BVf, VSI and snr

% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_dico       = 'files/dico_random_220000signals.mat';

permutations    = [1 2 3]; % by default: [BVf VSI snr] 

signal_number   = 200000;
test_number     = 2000; 

snr_max         = 100;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(file_dico)

K           = 15;
Lt          = 1;

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

Ytrain_noisy = Ytrain_noisy(:,permutations);
[theta, r]  = EstimateInverseFunction(Ytrain_noisy(:,1:2), Xtrain_noisy, K, Lt, 30, 0);


% Estimate parameters
Ytest_noisy = [Ytest zeros(length(Ytest), 1)];
Xtest_noisy = Xtest;
for i = 1:length(Xtest)
    Ytest_noisy(i,3)    = rand* snr_max;
    Xtest_noisy(i,:)    = AddNoise(Xtest(i,:), Ytest_noisy(i,3), 0);
end

Ytest_noisy = Ytest_noisy(:,permutations);
Ypredict 	= EstimateParametersFromModel(Xtest_noisy, theta);


%% Prediction analysis (correlation)

C = corr(Ypredict, Ytest_noisy);

figure
subplot(231);   plot(Ytest_noisy(:,1), Ypredict(:,1), '.')
title('Corr = ' +  string(C(1,1))); xlabel('Real values'); ylabel('Prediction')
subplot(232);   plot(Ytest_noisy(:,2), Ypredict(:,2), '.')
title('Corr = ' +  string(C(2,2)))
subplot(233);   plot(Ytest_noisy(:,3), Ypredict(:,3), '.')
title('Corr = ' +  string(C(3,3)))

YpredictNorm = Ypredict; % YpredictNorm is Ypredict but we put back values in their interval
YpredictNorm(Ypredict(:,1) > 0.25,:)    = NaN;
YpredictNorm(Ypredict(:,1) < 0,:)       = NaN;
YpredictNorm(Ypredict(:,2) > 1e-5,:)    = NaN;
YpredictNorm(Ypredict(:,2) < 1e-6,:)    = NaN;

CNorm = corr(YpredictNorm, Ytest_noisy, 'rows','pairwise');

subplot(234);   plot(Ytest_noisy(:,1), YpredictNorm(:,1), '.')
title('Corr = ' +  string(CNorm(1,1)))
subplot(235);   plot(Ytest_noisy(:,2), YpredictNorm(:,2), '.')
title('Corr = ' +  string(CNorm(2,2)))
subplot(236);   plot(Ytest_noisy(:,3), YpredictNorm(:,3), '.')
title('Corr = ' +  string(CNorm(3,3)))
% p = polyfit(Ytest_noisy(:,3), Ypredict(:,3), 1);
% hold on; plot(Ytest_noisy(:,3), p(1) * Ytest_noisy(:,3) + p(2), 'g.');
% 


