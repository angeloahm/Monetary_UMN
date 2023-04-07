function u = utility(c)

% Define utility function
if c<0
    u=-1e7;
else
    u = c.^(1-sigma) / (1-sigma);
end