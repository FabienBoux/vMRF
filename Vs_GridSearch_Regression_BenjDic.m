
addpath(genpath(fullfile(pwd, 'functions')))


%% Regression approach Vs Grid search approach

% Load data
file_dico_grid        = 'files/dicos_benj/dico_ratiopostpre_benj.mat';
file_dico_regression  = 'files/dicos_benj/dico_ratiopostpre_benj.mat';
file_dico_test        = 'files/dicos_benj/dico_all_S.mat';

lw      = 1;
coord   = [1 2];


%% Learning step

% Prepare dico reg
load(file_dico_regression)
Xreg    = abs(X);
Yreg    = Y(:,coord);
[Yreg, rmved] = rmmissing(Yreg);
Xreg    = Xreg(~rmved,:);

% Compute regression prediction
cstr.Sigma      = 'd';
[theta, ~]   	= EstimateInverseFunction(Yreg, Xreg, 50, lw, 200, cstr);


%% Test step

% Prepare dico grid
load(file_dico_grid)
Xgrid  	= X;
Ygrid  	= Y(:,coord);

% Prepare dico test
load(file_dico_test)
Xtest  	= X;
Ytest  	= Y;
vect    = (Ytest(:,1) < 10 ) & (Ytest(:,1) > 1 ) & (Ytest(:,2) > 1 ) & (Ytest(:,2) < 10 );
Ytest   = Ytest(vect,:);
Xtest   = Xtest(vect,:);

% Prepare reference data
Yf      = Yf(vect,:);
clear X Y

% Compute regression prediction
Ypredict_reg	= EstimateParametersFromModel(Xtest, theta);

% Compute reference predictions
Ypredict_grid 	= EstimateParametersFromGrid(Xtest, Xgrid, Ygrid);

Ypredict_grid(:,1) = 100* Ypredict_grid(:,1);
Ypredict_reg(:,1)  = 100* Ypredict_reg(:,1);
% Ypredict_grid(:,2) = 1e6* Ypredict_grid(:,2);
% Ypredict_reg(:,2)  = 1e6* Ypredict_reg(:,2);

% Remove absurd data
Ypredict_reg(Ypredict_reg(:,1) > 10) = nan;
Ypredict_reg(Ypredict_reg(:,1) <  1) = nan;
Ypredict_reg(Ypredict_reg(:,2) > 1e-5, 2) = nan;
Ypredict_reg(Ypredict_reg(:,2) <  1e-6, 2) = nan;

% Evaluation predictions
tmp = abs(Ytest - Ypredict_grid(:,1:2));
rslt(1,:) = mean(rmmissing(tmp));
tmp = abs(Ytest - Ypredict_reg(:,1:2));
rslt(2,:) = mean(rmmissing(tmp));
tmp = abs(Ytest - Yf);
rslt(3,:) = mean(rmmissing(tmp));


% Display
% figure
% subplot(131); plot(Ytest(:,1), Ypredict_reg(:,1),  '.')
% hold on; plot([0 10], [0 10], 'r')
% title(num2str(rslt(1,1)))
% 
% subplot(132); plot(Ytest(:,1), Ypredict_grid(:,1), '.')
% hold on; plot([0 10], [0 10], 'r')
% title(num2str(rslt(2,1)))
% 
% subplot(133); plot(Ytest(:,1), Yf(:,1), '.')
% hold on; plot([0 10], [0 10], 'r')
% title(num2str(rslt(3,1)))

%% Plot

figure
subplot(231); densityplot(Ytest(:,1), Ypredict_reg(:,1), [400 200])
hold on; plot([0 10], [0 10], 'k', 'LineWidth', 1.5)
xlim([1 10]); ylim([1 10])
title('Regression vs Analyse classique')

subplot(232); densityplot(Ytest(:,1), Ypredict_grid(:,1), [400 200])
hold on; plot([0 10], [0 10], 'k', 'LineWidth', 1.5)
xlim([1 10]); ylim([1 10]);
title('Dan Ma vs Analyse classique')

subplot(233); densityplot(Ytest(:,1), Yf(:,1), [400 400])
hold on; plot([0 10], [0 10], 'k', 'LineWidth', 1.5)
title('Benjamin vs Analyse classique')
xlim([1 10]); ylim([1 10]);


subplot(234); densityplot(Ytest(:,2), Ypredict_reg(:,2), [200 200])
hold on; plot([0  10], [1e-6  1e-5], 'k', 'LineWidth', 1.5)
xlim([1 10]); ylim([1e-6 1e-5]);

subplot(235); densityplot(Ytest(:,2), Ypredict_grid(:,2), [200 800])
hold on; plot([0  10], [1e-6  1e-5], 'k', 'LineWidth', 1.5)
xlim([1 10]); ylim([1e-6 1e-5]);

subplot(236); densityplot(Ytest(:,2), Yf(:,2), [200 800])
hold on; plot([0 10], [0 10], 'k', 'LineWidth', 1.5)
xlim([1 10]); ylim([1 10]);



