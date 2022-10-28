function [x, k] = uo_nn_GM(x, f, g, epsG, kmax, ils, ialmax, kmaxBLS, epsal, c1, c2)
%
% Input parameters: x : initial point. f : function to be minimized. g :
% gradient of the function to be minimized. epsG : optimality tolerance. kmax :
% maximum number of iterations. ils : line search (1 if exact, 2 if uo_BLS, 3 if
% uo_BLSNW32) ialmax :  formula for the maximum step lenght (1 or 2). kmaxBLS :
% maximum number of iterations of the uo_BLSNW32. epsal : minimum progress in
% alpha, algorithm up_BLSNW32 c1,c2 : (WC) parameters.
%
% Output parameters: x : optimal point. k : number of iterations.
%

% Initialization:
    k = 0;

    while (norm(g(x)) > epsG) && (k < kmax)
        % direction
        d = -g(x);

        % Line search:
        if ils == 1
            % Exact line search:
            Q = hess(f, x); a = -(Q*x)'*d/(d'*Q*d);
        elseif ils == 2
            % uo_BLS:
            [a, ~] = uo_BLS(f, g, x, d, ialmax);
        elseif ils == 3
            % uo_BLSNW32:
            [a, ~] = uo_BLSNW32(f, g, x, d, alpham, c1, c2, kmaxBLS, eps);
        else
            error('Error: ils must be 1, 2 or 3');
        end

        % Update:
        x = x + a*d; k = k + 1;
    end

end
