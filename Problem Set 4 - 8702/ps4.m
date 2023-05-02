%% Problem Set 4 - 8702
% Prof. Javier Bianchi 
% Authors: Lucas Belmudes, Hasan Cetin, Jason Hall, Angelo Mendes

clear 
clc

%% Parameters 
% Model parameters
beta = 0.96;        % Discount factor
phi = 0.2;          % Default cost
gamma = 5;          % Risk aversion
pi = 0.5;           % Sunspot probability
y = 1;              % Agg. endowment

%R = 1.02;
R = 1/beta;         % Interest rate

% Define utility function
u = @(c) c.^(1-gamma)./(1-gamma);

% Numerical parameters 
print = 100;
eps = 1e-8;
max_iter = 1000;
n_b = 300;          %gridpoints
b_max = 0.99;

% Compute grid for debt
b_grid = linspace(0,b_max,n_b);

% Compute default value
Vd = (u(y)-phi)/(1-beta);

% Compute V_R when the gvt. cannot borrow anymore
V_r_low = u(y-b_grid) + beta/(1-beta) * u(y);

% Allocate arrays
V_r_high = zeros(n_b,1);
V0 = ones(n_b,1)*u(y);
RHS = zeros(n_b,1);
q0 = ones(n_b,1)*beta;
policy_idx = zeros(n_b,1);
policy_b = zeros(n_b,1);

%% Exercise (i)
% Find threshold between the safe zone and the crisis zone

G = @(b) u(y-b) + beta*(u(y)/(1-beta)) - Vd;
options = optimset('Display','off');
b_low = fsolve(G,-0.5, options);

fprintf('The threshold of debt between the safe and crisis zone is: %0.4f \n', b_low)

%% Exercise (ii)

for iter = 1:max_iter
    
    % Find regions according to the guess
    crisis_states = (V0>=Vd) .* (V_r_low < Vd);

    %Compute Bellman equation
    V = max(V0, Vd);
    for i_b = 1:n_b
        b = b_grid(i_b);
        for i_b_prime = 1:n_b
            b_prime = b_grid(i_b_prime);
            RHS(i_b_prime) = u(y - b + q0(i_b_prime)*b_prime) + ...                              % u(c) for each choice of b'
                             beta*( crisis_states(i_b)*(pi*Vd+(1-pi)*V_r_low(i_b_prime)) +...    % crisis zone: compute expectation
                                   (1-crisis_states(i_b))*V(i_b_prime) );                        % default or safe zone
        end
        [V_r_high(i_b), policy_idx(i_b)] = max(RHS);
        policy_b(i_b) = b_grid(policy_idx(i_b));
    end
    
    % Compute prices
    default_states = (V_r_high<Vd);                     %update default, safe, and crisis regions according to the updated V_r_high
    crisis_states = (V_r_high>=Vd) .* (V_r_low < Vd)';
    delta = default_states + pi*crisis_states;          %update probability of default, accounting for multiplicity region
    q = (1-delta) / R;                                  %update bond price schedule
    
    % Compute norm
    norm = max( max(abs(V_r_high(:)-V0(:))), max(abs(q(:)-q0(:)))  );

    %Check convergence
    if norm>eps
        q0 = q;
        V0 = V_r_high;
        if mod(iter,print)==0
            fprintf('Error: %0.8f        Iteration: %0.0i \n', norm, iter)
        end
    else
        fprintf('Solution has been found! \n')
        fprintf('Error: %0.8f        Iteration: %0.0i \n', norm, iter)
        break
    end
end

plot(b_grid, V_r_high, LineWidth=2, Color="#009900")
hold on 
plot(b_grid, V_r_low, LineWidth=2, Color="#EDB120")
plot(b_grid, Vd*ones(n_b,1), LineWidth=2, Color="#CC0000")
xline(b_low, LineWidth=2, LineStyle='--', Color="#C0C0C0")
b_high = max(b_grid(crisis_states>0));
%xline(b_high, LineWidth=2, LineStyle='--', Color="#C0C0C0")
legend('V^+_R', 'V^-_R','V_D')
%ylim([min(V_r_high(:)) max(V_r_high(:))])

ylim([-20 0])
%xlim([0 0.9])
hold off

fprintf('Interval for crisis region is: [%0.4f, %0.4f] \n', b_low, b_high)
