function [Kbest, Kopti] = FindOptimalK(Xtrain, Ytrain, Xtest, Ytest, K, Lw, maxiter, cstr, verb)

if ~exist('Kmax','var'),    K = 15:25; end
if ~exist('Lw','var'),      Lw = 0; end
if ~exist('maxiter','var'), maxiter = 200; end
if ~exist('cstr','var'),    cstr.Sigma  = 'd'; cstr.Gammat = 'd'; cstr.Gammaw = ''; end
if ~exist('verb','var'),    verb = 1; end


%% Compute the Kmax regression


Nrmse   = zeros(length(K), size(Ytrain,2));

parfor_progress(length(K));
parfor k = 1:length(K)
    
    [theta,~] 	= EstimateInverseFunction(Ytrain, Xtrain, K(k), Lw, maxiter, cstr, 0);
    Ypredict	= EstimateParametersFromModel(Xtest, theta, 0);

    for i = 1:size(Ytrain,2)
        Ypredict(Ypredict(:,i) > max(Ytrain(:,i)),i) = nan;
        Ypredict(Ypredict(:,i) < min(Ytrain(:,i)),i) = nan;
    end
    
    [~,Nrmse(k,:),~] = EvaluateEstimation(Ytest, Ypredict, Ytrain);
    
    parfor_progress;
end
parfor_progress(0);


%% Find the best K

[~,l]   = min(mean(Nrmse,2));
Kbest   = K(l);
Kopti   = Kbest;

% % Compute some model selection criteria
% %L = ; %maximum value of the likelihood function for the model
% M = size(Xtrain,2); N = size(Ytrain,2); %M is the X size and N is the Y size
% k = (1:Kmax)*(1 + M + N*M + M*(M+1)/2 + N); %number of estimated parameters in the model
% % This formula comes from the GLLiM paper
% n = size(Xtest,1); %number of observations / sample size
% 
% %R2 adjusted
% %TODO
% 
% %AIC / AICc
% aic         = 2 * (k - log(L));
% %AIC = N*log(RSS/N) + 2k
% [~,Kopti.aic] = min(mean(aic,2));
% 
% aicc        = aic + 2*k*(k+1) / (n-k-1);
% [~,Kopti.aicc] = min(mean(aicc,2));
% 
% %BIC
% bic         = k*log(n) - 2*log(L);
% [~,Kopti.bic] = min(mean(bic,2));

%Mallow's Cp
%TODO


if verb == 1
    figure
    
    plot(mean(Nrmse,2), 'o-', 'LineWidth', 1.5)
    hold on
    line('YData', ylim, 'XData', [Kbest Kbest], 'Color','r')
    %line('YData', ylim, 'XData', [Kopti Kopti], 'Color','g')
end


