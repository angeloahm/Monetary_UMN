% This code simulates the model to generate Tables 4-5 in AMR(2022)

d_sim_aux = zeros(T,1);
a_sim_aux = zeros(T,1);
b_sim_aux = zeros(T,1);
q_sim_aux = zeros(T,1);
z_sim_aux = zeros(T,1);
y_sim_aux = zeros(T,1);


[row, iz] = find(z_grid==1);
z_sim_aux(1) = z_grid(iz);
a_sim_aux(1) = a_grid(n_a);
q_sim_aux(1) = q(n_a, iz);

for t=1:T-1
    
    %Shock
    shock = unifrnd(0,1);
    
    % Map the random draw into the transition matrix
    cum_sum = cumsum(P(iz,:));
    iz = sum(shock >= cum_sum)+1; 
    
    ia = find(a_grid==a_sim_aux(t));
    
    z_sim_aux(t+1) = z_grid(iz);    
    a_sim_aux(t+1) = policy_a(ia, iz);
    d_sim_aux(t+1) = policy_d(ia, iz);
    b_sim_aux(t+1) = policy_b(ia, iz);
    q_sim_aux(t+1) = q(ia, iz);
    y_sim_aux(t+1) = z_sim_aux(t+1)*(1-gamma0*d_sim_aux(t+1)^gamma1);
    
end

z_sim = z_sim_aux(1001:T);
a_sim = a_sim_aux(1001:T);
d_sim = d_sim_aux(1001:T);
b_sim = b_sim_aux(1001:T);
q_sim = q_sim_aux(1001:T);
y_sim = y_sim_aux(1001:T);


%CONSOLIDATE INTO FILE:
%%%%%%%%%%%%%%%%%%%%%%%

a_data = a_sim_aux(T-50:T);
d_data = d_sim_aux(T-50:T);
y_data = y_sim_aux(T-50:T);

simulated_data = [a_data, d_data, y_data];


csvwrite('simulated_data.csv',simulated_data)