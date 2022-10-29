function [x, k] = uo_nn_GM(x, f, g, h, epsG, kmax, ils, ialmax, kmaxBLS, epsal, c1, c2)
%
% Input parameters:
%
% x : initial point.
% f : function to be minimized.
% g : gradient of the function to be minimized.
% h : Hessian of the function to be minimized.
% epsG : optimality tolerance.
% kmax : maximum number of iterations.
% ils : line search (1 if exact, 2 if uo_BLS, 3 if uo_BLSNW32)
% ialmax :  formula for the maximum step lenght (1 or 2).
% kmaxBLS : maximum number of iterations of the uo_BLSNW32.
% epsal : minimum progress in alpha, algorithm up_BLSNW32
% c1,c2 : (WC) parameters.
%
% Output parameters: x : optimal point. k : number of iterations.
%

% Initialization:
    k = 0;
    almin = 1e-5;
    rho = 0.5;
    iW = 1; % Wolfe condition

    while (norm(g(x)) > epsG) && (k < kmax)
        % direction
        d = -g(x);

        gd = g(x)'*d;
        fx = f(x);

        if k == 0
            almax = 1;
        elseif ialmax == 1
            % maximum step length
            almax = a*gd_prev/gd;
        else
            % maximum step length
            almax = 2*(fx - fx_prev)/gd;
        end

        % Line search:
        if ils == 1
            % Exact line search:
            a = -h(x)'*d/(d'*h(x)*d);
        elseif ils == 2
            % uo_BLS:
            [a, ~] = uo_BLS(x, d, f, g, almax, almin, rho, c1, c2, iW);
        elseif ils == 3
            % uo_BLSNW32:
            [a, ~] = uo_BLSNW32(f, g, x, d, almax, c1, c2, kmaxBLS, epsal);
        else
            error('Error: ils must be 1, 2 or 3');
        end

        gd_prev = gd;
        fx_prev = fx;

        % Update:
        x = x + a*d;
        k = k + 1;
    end

end
