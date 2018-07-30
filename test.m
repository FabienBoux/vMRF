
file = 'files/multidico/2018-07-25-09:35_32-samples_5-parameters_80000-signals.mat';
coord = [1 2];
K = 20; Lw = 0; cstr.Sigma = 'd';

load(file)
[theta, r]   	= EstimateInverseFunction(Y(:,coord), X, K, Lw, 100, cstr, 0);


%%
ind = randi([1 length(X)], 1,1);
y = X(ind,:)';
verb = 0;

%x_samples       = 0:0.0001:0.3;
x_samples       = 0:0.0001:0.03;
x_samples       = [x_samples; (1:9/(length(x_samples)-1):10)*1e-6];
[x_dens,psi]    = gllim_inverse_dens(y,theta,ones(size(y)),x_samples,0);
x_predict  		= gllim_inverse_map(y,theta,0);

clear tmp
figure
x_s = 0:0.000001:0.25;
subplot(121)
for k=1:20, tmp(k,:) = normpdf(x_s,psi.mu(1,k),psi.S(1,1,k)); end
plot(x_s, sum(tmp))
title([num2str(x_predict(1)) ' - ' num2str(Y(ind,1))])

% clear tmp
% x_s = (1:0.0001:10) *1e-6;
% subplot(122)
% for k=1:20, tmp(k,:) = psi.alpha(k) * normpdf(x_s,psi.mu(2,k),psi.S(2,2,k)); end
% plot(x_s, sum(tmp))
% title([num2str(x_predict(2)) ' - ' num2str(Y(ind,2))])


