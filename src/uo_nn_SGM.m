function [wo, fo, niter] = uo_nn_SGM(gL, Xtr, ytr, epsG, kmax, sg_al0, sg_be, sg_ga, sg_emax, sg_ebest, sg_seed)
% Stochastic Gradient Method for Unconstrained Optimization
%
% parameters:
%
% gL: gradient function
% Xtr: training data
% ytr: training labels
% epsG: stopping criterion
% kmax: maximum number of iterations
% sg_al0: initial step size
% sg_be: step size decay
% sg_ga: step size increase
% sg_emax: maximum number of consecutive steps without improvement
% sg_ebest: number of best steps to consider
% sg_seed: random seed
%
% returns:
%
% wo: optimal weights
% fo: optimal function value
% niter: number of iterations

end
