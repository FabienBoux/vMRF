

% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_dico       = 'files/dico_random_40000signals_MGEFIDSEprepost.mat';

signal_number   = 10000;

snr             = [100 75 50 40 30 25 20 17 14 11 8];

repetition      = 10;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(file_dico)

for rep = 1:repetition
    rand_perm   = randperm(length(X), signal_number);

    Xtrain      = abs(X(rand_perm, :));
    Ytrain      = Y(rand_perm, :);

    % Create model
    Xtrain_noisy = Xtrain;

    for i = 1:length(snr)
        Xtrain_noisy    = AddNoise(Xtrain, snr(i), 0);
        cstr.Sigma      = [];
        [theta, ~]      = EstimateInverseFunction(Ytrain, Xtrain_noisy, 15, 0, 30, cstr, 0);

        err(i,rep)    	= mean(diag(sum(theta.Sigma,3)));
    end
end

%%

figure
plot(snr, mean(err,2), 'o')