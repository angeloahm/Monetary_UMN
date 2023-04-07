% This code computes the VFI for the partial default model

% Define initial guess for V0
V = zeros(n_a, n_z);
V0 = ones(n_a, n_z);
for i_a=1:n_a
    for i_z=1:n_z
        z = z_grid(i_z);
        
        V0(i_a, i_z) = (z ^ (1-sigma) / (1-sigma)) / (1-beta);
    end
end

% Initial price guess
q0 = ones(n_a, n_z)*(1/R);
q = ones(n_a, n_z)*(1/R);

% Allocate policy functions
policy_a = ones(n_a, n_z);
policy_index = ones(n_a, n_z);
policy_d = zeros(n_a, n_z);
policy_b = zeros(n_a, n_z);

%Allocate arrays
V_aux = zeros(n_a,1);
RHS = zeros(n_a,n_z);

for iter=1:max_iter
    
    %Initial guess 
    %q = q0;
    
    for i_a=1:n_a %for each debt 
        a = a_grid(i_a);
        for i_z=1:n_z %for each shock
            z = z_grid(i_z);
            
            for i_aprime=1:n_a %for each possible a'
                a_prime = a_grid(i_aprime);
                
                %Compute default
                def = ( a*(1-R*kappa*q0(i_aprime, i_z))/(z*gamma0*gamma1) )^(1/(gamma1-1));
                
                % truncate default policy
                if def>1
                    d=1;
                else
                    if def<0
                        d=0;
                    else
                        d=def;
                    end
                end
                
                %compute consumption
                c = z*(1-gamma0*d^gamma1) - a*(1-d) + q0(i_aprime, i_z)*(a_prime-R*kappa*d*a);
                if c<0
                    u = -1e+7;
                else
                    u = (c ^(1-sigma)) / (1-sigma);
                end
            
                V_aux(i_aprime) = u + beta*dot(P(i_z,:), V0(i_aprime,:));
                
            end %a_prime
            
            %Optimal savings and value:
            [value, i_policy] = max(V_aux);
            a_policy = a_grid(i_policy);
            V(i_a, i_z) = value;
            policy_a(i_a, i_z) = a_policy;
            policy_index(i_a, i_z) = i_policy;

            %Use asset policy to update default policy
            def = ( a*(1-R*kappa*q0(i_policy, i_z))/(z*gamma0*gamma1) )^(1/(gamma1-1));
                
            % truncate default policy
            if def>1
                d=1;
            else
                if def<0
                    d=0;
                else
                    d=def;
                end
            end 
            
            % Optimal default and borrowing
            policy_d(i_a, i_z) = d;
            policy_b(i_a, i_z) = a_policy-R*kappa*d*a;
            
        end %z
    end %a
    
    
    for i_a=1:n_a
       a = a_grid(i_a);
       for i_z=1:n_z
          
           i_policy=policy_index(i_a, i_z);
           
           %Price update 
           RHS(i_a, i_z) = (1/R)*(1-policy_d(i_a, i_z)+ ...
                 R*kappa*policy_d(i_a, i_z).* q0(i_policy, i_z));
           
           %update bond prices
           q(i_a, i_z) = (1-zeta)*q0(i_a, i_z) + ...
               zeta*dot(P(i_z, :), RHS(i_a, :));
           
       end
    end
    
    %Check convergence
    norm_VF = max(abs(V(:)-V0(:)));
    norm_q = max(abs(q(:)-q0(:)));
    
    if norm_VF<tol && norm_q<tol
        disp('Equilibrium has been found!')
        disp(iter)
        break
    else
        if mod(iter,50)==0
            fprintf('Iter = %.5f \n',iter)
            fprintf('Error VF = %.8f \n',norm_VF)
            fprintf('Error  q = %.8f \n',norm_q)
        end
        V0=V;
        q0=q;
    end
    
end %iter




