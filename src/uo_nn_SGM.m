function [wo, k] = uo_nn_SGM(w, L, gL, Xtr, ytr, Xte, yte, sg_al0, sg_be, sg_ga, sg_emax, sg_ebest, sg_seed)
% Stochastic Gradient Method for Unconstrained Optimization
%
% parameters:
%
% w: initial w
% L: loss function
% gL: gradient of loss function
% Xtr: training data
% ytr: training labels
% Xte: test data
% yt: test labels
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
% k: number of iterations

% set random seed
rng(sg_seed);

% number of columns
p = size(Xtr, 2);

% size of mini-batch
m = floor(sg_ga * p);
k_sg_e = ceil(p/m);
k_sg_max = sg_emax * k_sg_e;
e = 0;
s = 0;

lte_best = Inf;
k = 0;

a_sg = 0.01*sg_al0;
k_sg = floor(sg_be * k_sg_max);

while e <= sg_emax && s < sg_ebest
    P = randperm(p);
    for i = 0:k_sg_e-1
        % mini-batch
        S = P(i*m+1:(i+1)*m);
        X = Xtr(:, S);
        y = ytr(S);

        % direction
        d = -gL(w, X, y);

        if k <= k_sg
            a = (1 - k/k_sg) * sg_al0 + k/k_sg * a_sg;
        else
            a = a_sg;
        end

        % update
        w = w + a * d;
        k = k + 1;
    end
    e = e + 1;
    lte = L(w, Xte, yte);

    if lte < lte_best
        lte_best = lte;
        % avg of each row in w
        wo = w;

        s = 0; % Break out of the loop
    else
        s = s + 1;
    end

end


end
