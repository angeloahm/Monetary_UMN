% Problem Set - Monetary Economics 8701
% Angelo Avelar Hermeto Mendes

% This code solves for the equilibrium of the simplified version of the 
% partial default model we saw in class

clear; clc;

tic
% Discretize shocks and set parameters
set_parameters

% Solve model by gridsearch 
solve_model

% Simulate the model 
simulate 

% Print results 
results

% Plot results
%plots


toc
