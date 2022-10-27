Q = [4 0; 0 1]; f = @(x) (1/2)*x'*Q*x; g = @(x) Q*x;
x = [1;2]; d = [-4;-2];
almax= 1.0; almin= 10^-6; rho = 0.5; c1 = 0.1; c2 = 0.5; iW=1;

[al,iWout] = uo_BLS(x,d,f,g,almax,almin,rho,c1,c2,iW)

% [start] Alg. BLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [al,iWout] = uo_BLS(x,d,f,g,almax,almin,rho,c1,c2,iW)
% iWout = 0: al does not satisfy any WC
% iWout = 1: al satisfies (WC1)
% iWout = 2: al satisfies WC
% iWout = 3: al satisfies SWC


end
% [end] Alg. BLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%