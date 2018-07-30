function [Yestim] = EstimateParametersFromModel(X, f_estim, verb)

narginchk(2, 3);
if nargin == 2, verb = 0; end


Yestim = gllim_inverse_map(X', f_estim, verb)';