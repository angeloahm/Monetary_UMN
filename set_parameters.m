
% Parameters
n_z = 5;
n_a = 50;
z_bar = 0;
rho = 0.875;
sigma_eta = 0.052;
sigma = 2;          %CRRA 
beta = 0.954;
R = 1/beta;
kappa = 0.7;
tol = 1e-5;
zeta = 0.1;
T = 30000;

max_iter = 2000;

% Default cost parameters 
gamma0 = 0.0476;
gamma1 = 2;

% Define the grid for z
[log_z_grid, P] = tauchen(n_z, 1, z_bar, rho, sigma_eta);
z_grid = exp(log_z_grid);

% Define the grid for a
a_grid = linspace(0, 0.5*exp(z_bar), n_a);








