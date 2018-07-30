
addpath(genpath(fullfile(pwd, 'functions')))


%% Verify that signals on grid are correctly estimate
% The error sum is null = it's ok = every signal on the grid is estimate by
% itself

% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_dico = 'files/04-Apr-2018_64-samples_4-parameters_185367-signals.mat';

dico_size = [50000:25000:150000];
signal_test_number = 40000;
snr_value = 10;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(file_dico)


% For each dictionary size
for iter = 1:length(dico_size)

    fprintf('Dictionary %d / %d \t %d signals \n', iter, length(dico_size), dico_size(iter))

    rand_perm   = randperm(length(X), dico_size(iter));

    Xtrain      = abs(X(rand_perm, :));
    Ytrain      = Y(rand_perm, :);

    Xtest       = abs(Xtrain(1:signal_test_number, :));
    Ytest       = Ytrain(1:signal_test_number, :);

    Xtest_noisy = AddNoise(Xtest , snr_value, 0);

    Ypredict  	= EstimateParametersFromGrid(Xtest_noisy, Xtrain, Ytrain);

    Error{iter} = abs(Ypredict - Ytest);
end


%% Plot the error map

for iter = 1:length(dico_size)
    errors = [Ytest, Error{iter}/max(Error{iter})]; 
    
    figure
    h = histogram(Ytest(:,1), 100);
    x = h.BinEdges;

    h = histogram(Ytest(:,2), 100);
    y = h.BinEdges;
    close

    error_map   = zeros(length(x), length(y));
    c           = error_map;

    for i = 1:length(errors)
        t = errors(i,:);

        [~, x_loc] = min(abs(x - t(1)));
        [~, y_loc] = min(abs(y - t(2)));

        error_map(x_loc, y_loc) = error_map(x_loc, y_loc) + t(3);
        c(x_loc, y_loc) = c(x_loc, y_loc) +1;
    end
    error_map = error_map ./ c;
    error_map(isnan(error_map)) = 0;

    [xx, yy] = meshgrid(x,y);

%     subplot(1,2,iter); surf(xx, yy, error_map)
    figure
    imagesc(x, y, error_map)
    title(string(dico_size(iter)))
end


    
