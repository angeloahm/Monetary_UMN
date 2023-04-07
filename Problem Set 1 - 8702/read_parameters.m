%% 1. Model Parameter Values

% MODEL PARAMETERS

sigma    = 2;
ita      = (1/0.83)-1; % elasticity parameter (elasticity is 0.83)
beta     = 0.906;       % discount factor
omega    = 0.3070;       % share for tradables
kappa    = 0.3235;
r        = 0.04;  R=1+r;


sep_prefs = 0;
if sep_prefs==1  % allows analytical solution of b*
    sigma =ita+1;  
end

NT_shock = 1;  % to use T-NT process for paper

str = sprintf('MODEL PARAMETERS: sigma %2.1f  elasticity %2.2f  beta %2.2f omega %2.2f kappa %2.2f  %r',sigma, 1/(1+ita),beta,omega,kappa,r); disp(str)

%%  2.  Numerical Parameters

if NT_shock ==1
    
    load proc_shock yT yN Prob 
    
    NSS =length(yT);
    
else
    
    yn          = 1   ;
    sigma_y     = 0.058; % stdev  yT
    rho         = 0.5;
    sigma_eps   = sqrt(sigma_y^2*(1-rho^2)); % stdev of innovation to yT
    
    NSS          = 9;
%     [Z,Prob]    = tauchenhussey(NSS,0,rho,sigma_eps ,sigma_eps );
     m=2;   
    [Z,Zprob] = tauchen(NSS,0,rho,sigma_eps,m);
    yT          = exp(Z);
    yN          = repmat(yn,NSS,1);
  
end

NB          = 100;
 

bmin_NDL    = -kappa*(1+min(yT(1)))+0.01;  % maximum debt consistent with feasibility
bmin        = -1.02;
bmax        = -0.2000;

B           = linspace(bmin,bmax,NB)';


SR          = ones(NB,NSS)*(1+r);
STAU        = ones(NB, NSS)*(1+tax);
YT          = repmat(yT,1,NB )';
YN          = repmat(yN,1,NB )';
KAPPAS      = kappa*ones(NB,NSS);
b           = repmat(B,1,NSS);
 

% ininitialize
bp          = b;
c           = STAU.*b.*SR  +YT  -bp +(1-STAU).*b;
 


price       = (1-omega)/omega*c.^(1+ita);
EE          = zeros(NB,NSS);
bmax_collat = -kappa*(price.*YN +YT);
cbind       = STAU.*b.*SR  +YT  -bmax_collat+(1-STAU).*b;
emu         = zeros(NB,NSS)*nan;
Ibind       = ones(NB,NSS);
 
 
uptd        = 0.2;   % Change to 0.1 if policy function fluctuate
uptdSP      = 0.02;

outfreq     = 50;    % Display frequency
iter_tol    = 500;   % Iteration tolerance
Tol_eulerDE = 1.0e-7;% when to look for unconstrained solution
Tol_eulerSP = 1.0e-7;
tol         = 1.0e-5; %Numerical tolerance


T   = 10000;  % Number of simulations
cut = 0;  % Initial condition dependence


options = optimset('Display','off');
str     = sprintf('tol %10.2e  bmin %2.2f bmax %2.2f  NB %3i NS %2i',tol,bmin,bmax,NB,NSS); disp(str)

disp(' ')