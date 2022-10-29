function [x, k] = uo_nn_QNM(x, f, g, epsG, kmax, ils, ialmax, kmaxBLS, epsal, c1, c2)
    % Quasi-Newton method / BFGS Update

    I = eye(length(x));
    H = I;
    k = 0;
    iW = 1;
    almin = 1e-5;
    rhoBLS = 0.5;

    while norm(g(x)) > epsG && k < kmax
        d = -H*g(x);

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

        gd_prev = gd;
        fx_prev = fx;

        if ils == 1
            % Exact line search:
            a = -h(x)'*d/(d'*h(x)*d);
        elseif ils == 2
            % uo_BLS:
            [a, ~] = uo_BLS(x, d, f, g, almax, almin, rhoBLS, c1, c2, iW);
        elseif ils == 3
            % uo_BLSNW32:
            [a, ~] = uo_BLSNW32(f, g, x, d, almax, c1, c2, kmaxBLS, epsal);
        else
            error('Error: ils must be 1, 2 or 3');
        end

        x_prev = x;
        x = x + a*d;
        s = x - x_prev;

        y = g(x) - g(x_prev);
        rho = 1/(y'*s);

        H = (I - rho*s*y')*H*(I - rho*y*s') + rho*s*s';
        k = k + 1;
    end
end
