%% Exercise 2: Fixed Tax Experiment

clc;
close all

addpath('auxiliar')

%% Initialize code

tax=0;

% Read model/numerical parameters
read_parameters

% Solve for the descentralized eq.
equilibrium

% Solve for the SPP
social_planner

% Compare welfare 
compare_welfare

V_de = V;

%Simulate model
simulation

sd_CASP = std(CAYSP_SIM(:));                                %SD of the current account
crises = (CAYSP_SIM>sd_CASP) .* bindSP_SIM;                   %Define the crise event
freq_crisesSP = sum(crises(:))/T;                           %Compute freq. 
avg_CASP = sum(CAYSP_SIM(:) .* crises(:)) / sum(crises(:));     %Compute avg. CA during crises

%Compute avg. tau
avg_opt_tax = sum(Tau_SIM(:)) / T;


%%% GRID FOR TAU 
n_tau=50;
tau_grid = linspace(0, 0.1, n_tau);
%Egamma_vec = zeros(1,n_tau);
max_gamma_vec = zeros(1,n_tau);
gamma_matrix = zeros(NB,n_tau);
b_where = zeros(1,n_tau);
freq_crises_vec = zeros(1,n_tau);
avg_CA_vec = zeros(1,n_tau);


for t=1:length(tau_grid)
    tax = tau_grid(t);
    fprintf('=========================\n')
    fprintf('Solving for tau=%2f',tax);
    fprintf('\n')
    fprintf('=========================\nd')
    
    % Solve for the equilibrium when tau>0
    read_parameters
    equilibrium
    compare_welfare
    

    % Compute welfare gains 
    gamma = (V ./ V_de).^(1/(1-sigma))-1;

    %Compute stationary distribution
    simulation
    %distribution
    
    % Compute welfare gap and find where 
    gamma_matrix(:,t) = gamma(:,6);
    %E_gamma = sum(Pi(:) .* gamma(:));
    %Egamma_vec(t) = E_gamma;
    [max_gamma, idx] = max(gamma(:,6));
    max_gamma_vec(t) = max_gamma;
    b_where(t) = B(idx);
    
    % Compute crises statistics
    sd_CA = std(CAY_SIM(:));                                    %SD of the current account
    crises = (CAY_SIM>sd_CA) .* bindSIM;                        %Define the crise event
    freq_crises = sum(crises(:))/T;                             %Compute freq. 
    freq_crises_vec(t) = freq_crises;
    
    avg_CA = sum(CAY_SIM(:) .* crises(:)) / sum(crises(:));     %Compute avg. CA during crises
    avg_CA_vec(t) = avg_CA;

end

figure(1)
plot(tau_grid, max_gamma_vec, LineWidth=2)
title('Welfare gains for given (\tau)')
xlabel('\tau')
ylabel('Welfare gain')
saveas(gcf,'welfare_gain.png')


figure(2)
plot(tau_grid, b_where, LineWidth=2)
title('Debt that maximizes welfare gains')
xlabel('\tau')
ylabel('B')
saveas(gcf,'bond_holdings.png')


figure(3)
plot(tau_grid, freq_crises_vec*100, LineWidth=2)
title('Frequency of crises')
xlabel('\tau')
ylabel('Freq. of crises')
saveas(gcf,'frequency.png')

figure(4)
plot(tau_grid, avg_CA_vec, LineWidth=2)
title('Avg. CA in crisis')
xlabel('\tau')
ylabel('Current account')
saveas(gcf,'CA.png')

% figure(3)
% plot(B, gamma_matrix, LineWidth=2)
% title('Welfare gains for given (\tau)')
% xlabel('B')
% ylabel('Welfare gain')














