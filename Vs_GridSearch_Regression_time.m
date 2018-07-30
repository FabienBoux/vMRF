addpath(genpath(fullfile(pwd, 'functions')))

%% Regression approach Vs Grid search approach

% Load data
file_dico_grid        = 'files/dico_grid_216931signals.mat';
file_dico_regression  = 'files/dico_random_220000signals.mat';

dic_min     = 1;
dic_max     = 10;
nb_signal   = 10000;
repetition  = 10;
coord       = [1 2];

lw          = 0;
cstr.Sigma  = 'd';

load(file_dico_regression)
Xreg        = abs(X);
Yreg        = Y(:,coord);

load(file_dico_grid)
Xgrid       = abs(X);
Ygrid       = Y(:,coord);
clear X Y

load(file_dico_regression)


% Grid search
for rep = 1:repetition
    fprintf(['Repetition: \t ' num2str(rep) ' / ' num2str(repetition) '\n'])
    
    % Select test signals
    rand_perm   = randperm(length(Xreg), nb_signal);
    Xtest       = Xreg(rand_perm,:);
    Ytest       = Yreg(rand_perm,:);
    
    
    
    for dic = dic_min:dic_max
        fprintf(['\t Grid \t\t - Iteration: \t ' num2str(dic) ' / ' num2str(dic_max) '\n'])

        vect            = zeros(1, length(unique(Ygrid(:,1))));
        vect(1:dic:end) = 1;

        zero            = zeros(size(vect));
        vect            = repmat([vect repmat(zero, 1, dic-1)], 1, ceil(length(Ygrid) / length([vect repmat(zero, 1, dic-1)])));
        vect            = logical(vect);

        Xgrid_eff       = Xgrid(vect,:);
        Ygrid_eff       = Ygrid(vect,:);
        
        dico_len(dic)   = length(Ygrid_eff);
        
        tic;
        Ypredict_grid   = EstimateParametersFromGrid(Xtest, Xgrid_eff, Ygrid_eff);
        grid_time(dic, rep) = toc;
    end
    

    % Regression
    for dic = dic_min:dic_max
        fprintf(['\t Regression \t - Iteration: \t ' num2str(dic) ' / ' num2str(dic_max) '\n'])
        
        rand_perm       = randperm(length(Xreg), dico_len(dic));
        Xregression_eff = Xreg(rand_perm,:);
        Yregression_eff = Yreg(rand_perm,:);        
        
        tic
        [theta, ~]    	= EstimateInverseFunction(Yregression_eff, Xregression_eff, 20, lw, 100, cstr);
        learning_time(dic, rep) = toc;
        
        tic
        Ypredict_reg    = EstimateParametersFromModel(Xtest, theta);
        
        for i = 1:size(Ygrid,2)
            Ypredict_reg(Ypredict_reg(:,i) > max(Ygrid(:,i))) = nan;
            Ypredict_reg(Ypredict_reg(:,i) < min(Ygrid(:,i))) = nan;
        end
        reg_time(dic, rep) = toc;
    end
end


%% Display
figure
subplot(121); errorbar(dico_len, mean(grid_time'), std(grid_time'), 'o-', 'LineWidth', 1.5)
hold on
errorbar(dico_len, mean(reg_time'), std(reg_time'), 'o-', 'LineWidth', 1.5)
ylabel('Time (in s)'); xlabel('Dictionary length')
legend({'Grid approach', 'Regression approach'})
title(['Evaluation of ' (num2str(nb_signal)) ' ' num2str(size(Xgrid,2)) 'sample-signals'  ' observations on ' num2str(repetition) ' repetitions'])
set(gca,'XScale','log');set(gca,'YScale','log');

subplot(122); errorbar(dico_len, mean(learning_time'), std(learning_time'), 'o-', 'LineWidth', 1.5)
ylabel('Time (in s)'); xlabel('Dictionary length')
legend({'Learning step'})
title(['Evaluation of ' (num2str(nb_signal)) ' ' num2str(size(Xgrid,2)) 'sample-signals'  ' observations on ' num2str(repetition) ' repetitions'])
set(gca,'XScale','log'); set(gca,'YScale','log');


figure
modelfun    = @(a,x)a(1) + a(2) * x;
mdl         = fitnlm(dico_len, mean(learning_time'), modelfun, [1 1e-10]);
a           = mdl.Coefficients.Estimate;

v = min(dico_len):100:max(dico_len);
hold on
plot(dico_len, mean(grid_time'), 'ko-')
plot(dico_len, mean(reg_time'), 'kx-')
plot(v, modelfun(a,v), '-', 'LineWidth', 2)
plot(v, mean(mean(reg_time')) * ones(size(v)), '-', 'LineWidth', 2)
set(gca,'XScale','log'); set(gca,'YScale','log');
legend({'Grid approach', 'Learning approach', 'Grid linear regression', 'Learning linear regression'})
title(['Evaluation of ' (num2str(nb_signal)) ' ' num2str(size(Xgrid,2)) 'sample-signals'  ' observations on ' num2str(repetition) ' repetitions'])
