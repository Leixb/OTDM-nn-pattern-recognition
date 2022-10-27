% [start] Alg. BLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [al, iWout] = uo_BLS(x, d, f, g, almax, almin, rho, c1, c2, iW)
  % iWout = 0: al does not satisfy any WC
  % iWout = 1: al satisfies (WC1)
  % iWout = 2: al satisfies WC
  % iWout = 3: al satisfies SWC

  % input validation
  if ~ismember(iW, [1, 2])
    error('iW must be 1 (WC) or 2 (SWC)')
  end

  % Initialize
  al = almax;
  iWout = 0;

  while (al > almin)
    if (f(x + al*d) <= f(x) + c1*al*g(x)'*d) % (WC1)
      % (WC1) is satisfied, let's check (WC2 or SWC)
      if (iW == 1) % check WC2
        if (g(x + al*d)'*d >= c2*g(x)'*d)
          iWout = 2;
        end
        iWout = 1;
        return;
      else % check SWC
        if (abs(g(x + al*d)'*d) <= c2*abs(g(x)'*d))
          iWout = 3;
          return;
        end
      end
    end
    al = rho*al; % Adjust step size
  end
end
