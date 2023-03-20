% This code prints the results from Tables 4-5 in AMR after simulating

% Default results
freq_default = sum(d_sim>0);
cond_avg_default = mean(d_sim(d_sim>0));
cond_sd_default = std(d_sim(d_sim>0));

fprintf(' PARTIAL DEFAULT \n')
fprintf('Frequency of Default:      %.5f \n',100*freq_default/T)
fprintf('Avg. Default | Default>0:  %.5f \n',100*cond_avg_default)
fprintf('SD. Default | Default>0:   %.5f \n',100*cond_sd_default)

% Debt to output
avg_debt_output = mean(b_sim ./ y_sim);
sd_debt_output = std(b_sim ./ y_sim);
fprintf(' DEBT TO OUTPUT \n')
fprintf('Avg. Debt to Output:      %.5f \n',avg_debt_output*100)
fprintf('SD. Debt to Output:       %.5f \n',sd_debt_output*100)

% Debt service
avg_service_output = mean((b_sim ./ q_sim) ./ y_sim);
sd_service_output = std((b_sim ./ q_sim) ./ y_sim);
fprintf(' SERVICE TO OUTPUT \n')
fprintf('Avg. Service to Output:      %.5f \n',avg_service_output*100)
fprintf('SD. Service to Output:       %.5f \n',sd_service_output*100)

% Debt due
fprintf(' DEBT DUE \n')
avg_due_output = mean(a_sim ./ y_sim);
fprintf('Avg. Debt due to output:    %.5f \n',100*avg_due_output)


% Spread 
fprintf(' SPREAD \n')
spread = 100*(1./q_sim(q_sim>0) - R);
avg_spread = mean(spread);
sd_spread = std(spread);
fprintf('Avg. Spread:        %.5f \n',avg_spread)
fprintf('SD. Spread:         %.5f \n',sd_spread)
fprintf('Corr(y, spread):    %.5f \n',corr(spread, y_sim))
fprintf('Corr(debt, spread): %.5f \n',corr(spread, b_sim))

% Output 
fprintf(' OUTPUT \n')
mdl = arima(1,0,0);
mdl.Constant = 1-rho;

estimate_output = estimate(mdl, y_sim);
persistence = estimate_output.AR{1};
sd_output = sqrt(estimate_output.Variance);
fprintf('Persistence:   %.5f \n',persistence)
fprintf('SD. Output:    %.5f \n',100*sd_output)


estimate_shock = estimate(mdl, z_sim);
estimated_rho = estimate_shock.AR{1};
fprintf('Persistence (rho):   %.5f \n',estimated_rho)


% DEFAULTED COUPONS
fprintf(' DEFAULTED COUPONS \n')
defaulted = d_sim .* a_sim;
avg_defaulted = mean(defaulted(d_sim>0));
std_defaulted = std(defaulted(d_sim>0));
fprintf('Avg Defaulted | def>0:   %.5f \n',100*avg_defaulted)
fprintf('SD. Defaulted | def>0:   %.5f \n',100*std_defaulted)

% CONSUMPTION
fprintf(' CONSUMPTION \n')
c_sim = y_sim - a_sim.*(1-d_sim) + q_sim.*b_sim;
fprintf('SD. c (relative to output):   %.5f \n', std(c_sim)/std(y_sim))



