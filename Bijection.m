
file = 'dictionaries/multidico/18-08-23/2018-08-22-18:30_32-samples_4-parameters_6750-signals.mat';

load(file)

% 100 signals is enought
nb_signals = 1000;
ind     = randperm(size(X,1), nb_signals);
X       = X(ind,:);
Y       = Y(ind,:);

% Ratio dico
X       = X(:,end/2+1:end) ./ X(:,1:end/2);

% Dot product on normalized signals (correlation is an alternative)
Xnorm   = (1 ./ sum(X.^2, 2) .^0.5) * ones(1, size(X, 2)) .* X;
dx      = Xnorm*Xnorm'; %dx      = corr(X');
for i=1:size(X,1), for j=1:size(X,1), Rmse(i,j,:) = EvaluateEstimation(X(i,:), X(j,:), X); end; end
dx = mean(Rmse,3);

Ynorm   = Y; 
for c = 1:size(Y,2), Ynorm(:,c) = (Y(:,c) - min(Y(:,c))) ./ max(Y(:,c)); end
Ynorm   = (1 ./ sum(Ynorm.^2, 2) .^0.5) * ones(1, size(Ynorm, 2)) .* Ynorm;
dy      = Ynorm*Ynorm'; 

for i=1:size(Y,1), for j=1:size(Y,1), [~,Nrmse(i,j,:),~] = EvaluateEstimation(Y(i,:), Y(j,:), Y); end; end
dy      = mean(Nrmse,3);

figure
plot(reshape(dx,1,[]),reshape(dy,1,[]),'.')
xlabel('dX'); ylabel('dY')