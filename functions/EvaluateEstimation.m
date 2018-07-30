function [Rmse, Nrmse, Mae] = EvaluateEstimation(Ytrue, Yestim, Ytrain)

N = size(Ytrue,2);

Rmse = zeros(1, N);
Nrmse = Rmse;

for i = 1:N
    Rmse(i)  = nanmean( ( Ytrue(:,i) - Yestim(:,i) ).^2 )^.5;
    Nrmse(i) = Rmse(i) / nanmean( ( Ytrue(:,i) - nanmean(Ytrain(:,i)) ).^2 ).^0.5;
    Mae(i)   = nanmean( abs(Ytrue(:,i) - Yestim(:,i)) );
end
