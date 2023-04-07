% Function that uses Tauchen's method to discretize an AR(1) process.
function [grid, P] = tauchen(n, m, mu, rho, sigma)

% lb: grid's lower bound
% ub: grid's upper bound
% n: number of gridpoints
% m: scaling parameter
% mu: mean
% rho: AR(1) coefficient (autocorrelation)
% sigma: SD of AR process

% Defining variables
ub = mu+m*sigma/((1-rho^2)^(1/2)); % define upper bound
lb = mu-m*sigma/((1-rho^2)^(1/2)); % define lower bound

grid = linspace(lb, ub, n);     % equidistant points

delta = (ub-lb)/(n-1);          % distance between two points

P = zeros(n,n);

for i=1:n
    for j=2:n-1
        x_upper = (grid(j)+delta/2-(1-rho)*mu-rho*grid(i))/sigma;
        x_lower = (grid(j)-delta/2-(1-rho)*mu-rho*grid(i))/sigma;
        P(i,j) = normcdf(x_upper)-normcdf(x_lower);
    end
    x_1 =  (grid(1)-(1-rho)*mu-rho*grid(i)+delta/2)/sigma;
    P(i,1) = normcdf(x_1);
    
    x_n = (grid(n)-(1-rho)*mu-rho*grid(i)-delta/2)/sigma;
    P(i,n) = 1-normcdf(x_n);
end

