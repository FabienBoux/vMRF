function [Estimation, Parameters, Errors] = EstimateParametersBothMethod(Xtrain, Ytrain, Xgrid, Ygrid, Xtest, Ytest)

% Parameters
Parameters.maxiter     = 200;
Parameters.Lw          = 0;

Parameters.cstr.Sigma  = 'd';
Parameters.cstr.Gammat = 'd';
Parameters.cstr.Gammaw = '';

if ~exist('Ytest','var')
    Parameters.K = 20;
else
    % TODO: Implement automatic search of K (and maybe Lw ?)
    K       = 15:40;
    Parameters.K = FindOptimalK(Xtrain, Ytrain, Xtest, Ytest, K, Parameters.Lw, Parameters.maxiter, Parameters.cstr, 0);
end


% Estimations
Parameters.theta = EstimateInverseFunction(Ytrain, Xtrain, Parameters.K, Parameters.Lw, Parameters.maxiter, Parameters.cstr, 0);
Estimation.Regression.Y = EstimateParametersFromModel(Xtest, Parameters.theta, 0);

for i = 1:size(Ytrain,2)
    Estimation.Regression.Y(Estimation.Regression.Y(:,i) > max(Ygrid(:,i)),i) = nan;
    Estimation.Regression.Y(Estimation.Regression.Y(:,i) < min(Ygrid(:,i)),i) = nan;
end

Estimation.GridSearch.Y = EstimateParametersFromGrid(Xtest, Xgrid, Ygrid);
Estimation.Real.Y       = Ytest;

% Evaluation
[Errors.GridSearch.Rmse, Errors.GridSearch.Nrmse, Errors.GridSearch.Mae, Errors.GridSearch.Nmae] = EvaluateEstimation(Ytest, Estimation.GridSearch.Y, Ygrid);
[Errors.Regression.Rmse, Errors.Regression.Nrmse, Errors.Regression.Mae, Errors.GridSearch.Nmae] = EvaluateEstimation(Ytest, Estimation.Regression.Y, Ytrain);