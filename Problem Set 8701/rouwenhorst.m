function [grid, P] = rouwenhorst(n, mu, rho, sigma)

% lb: grid's lower bound
% ub: grid's upper bound
% n: number of gridpoints
% mu: drift
% rho: AR(1) coefficient (autocorrelation)
% sigma: SD of AR process

% Defining variables
var_theta = (sigma^2) / (1-rho^2);
sigma_theta = sqrt(var_theta);
ub = sigma_theta*sqrt(n-1);         %defining upper bound
lb = -ub;                           %defining lower bound

grid = linspace(lb, ub, n);         % equidistant points

p = (1+rho)/2;                      %defining p

% Create a cell array to store matrices with different dimensions
C = cell([1,n]);

C{1} = 1;
%P_2 = [p, 1-p; 1-p, p]
    
for i=2:n
    zero_vec = zeros(i-1,1);    %create zeros vector
    % Iterate prices
    C{i} = p*[C{i-1}, zero_vec; zero_vec', 0]+(1-p)*[zero_vec, C{i-1}; 0, zero_vec']+(1-p)*[zero_vec',0;C{i-1}, zero_vec]+p*[0,zero_vec';zero_vec, C{i-1}];

end

P = C{n};

for i=1:n
    P(i,:) = P(i,:)/sum(P(i,:));
end
    
