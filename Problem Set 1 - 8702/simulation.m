%% 7. Simulation

% load shocks_initials_aer.mat  b0_DE      b0_SP      state_sim
% S_index=state_sim;
disp('simulating the model')
% cut=0;
b0_DE     = -0.9443;
b0_SP     = -0.9443;
[S_index] = markov(Prob,T+cut+1,1,1:length(Prob));

% SHOCKS
ySIM     = YT(1,S_index)';
yNSIM    = YN(1,S_index)';
RSIM     = SR(1,S_index)';
KAPPASIM = KAPPAS(1,S_index)';

% INITIALIZING MATRIXE
bpSIM       = zeros(T,1);
WfSIM       = zeros(T,1);
bindSIM     = zeros(T,1);
cSIM        = zeros(T,1);
pSIM        = zeros(T,1);
CA_SIM      = zeros(T,1);
bpSP_SIM    = zeros(T,1);
bindSP_SIM  = zeros(T,1);
cSP_SIM     = zeros(T,1);
pSP_SIM     = zeros(T,1);
CASP_SIM    = zeros(T,1);
Tau_SIM     = zeros(T,1);
Tregion_SIM = zeros(T,1);

%Initial Conditions
bpSIM(1)    = b0_DE;
bpSP_SIM(1) = b0_SP;

for i=1:T
    
    bptemp     = interp1(B,bpDE(:,S_index(i)),bpSIM(i),'linear','extrap');
    cSIM(i)    = bpSIM(i)*RSIM(i)+ySIM(i)-bptemp;
    
    WfSIM(i)   = interp1(B,Wf(:,S_index(i)),bpSIM(i),'linear','extrap');
    
    bptempSP   = interp1(B,bpSP(:,S_index(i)),bpSP_SIM(i),'linear','extrap');
    cSP_SIM(i) = bpSP_SIM(i)*RSIM(i)+ySIM(i)-bptempSP;
    Tau_SIM(i) = interp1(B,tau(:,S_index(i)),bpSP_SIM(i),'linear','extrap');
    if i<T
        bpSIM(i+1)    = bptemp;
        bpSP_SIM(i+1) = bptempSP;
    end
end

fprintf('welfare gain is: ');
fprintf('%.6f \n ',  mean(WfSIM)*100);
disp(' ')

pSIM               = (1-omega)/omega*(cSIM./yNSIM).^(1+ita);
pSP_SIM            = (1-omega)/omega*(cSP_SIM./yNSIM).^(1+ita);

bplim_SIM0         = -  KAPPASIM.*(pSIM.*yNSIM+ySIM);
bplimSP_SIM0       = -  KAPPASIM.*(pSP_SIM.*yNSIM+ySIM);

bplim_SIM          = bplim_SIM0;
bplimSP_SIM        = bplimSP_SIM0;

bplim_SIM(2:end)   = bplim_SIM0(1:end-1);
bplimSP_SIM(2:end) = bplimSP_SIM0(1:end-1);

%%%%% IMPORTANT TO USE TOLERANCE NOT <= %%%%%%%%%%%%
bindSIM((bpSIM-bplim_SIM)<=1e-5) = 1;
bindSP_SIM((bpSP_SIM-bplimSP_SIM)<=1e-5) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bindSIM(1)    = 0; % Not binding the first perdiod
bindSP_SIM(1) = 0; % Not binding the first perdiod

bpSIM(bindSIM==1)       = bplim_SIM(bindSIM==1);
bpSP_SIM(bindSP_SIM==1) = bplimSP_SIM(bindSP_SIM==1);

Tregion_SIM(Tau_SIM>0) = 1;

CA_SIM(1:end-1)        = bpSIM(2:end)-bpSIM(1:end-1);  % >0 capital outflow
CASP_SIM(1:end-1)      = bpSP_SIM(2:end)-bpSP_SIM(1:end-1);

CA_SIM(T)              = interp1(B,bpDE(:,S_index(T)),bpSIM(T),'linear','extrap')-bpSIM(T);
CASP_SIM(T)            = interp1(B,bpSP(:,S_index(T)),bpSP_SIM(T),'linear','extrap')-bpSP_SIM(T);

Y_SIM                  = pSIM.*yNSIM+ySIM;   % Income
YSP_SIM                = pSP_SIM.*yNSIM+ySIM;

CAY_SIM                = CA_SIM./Y_SIM;  % Current Account
CAYSP_SIM              = CASP_SIM./YSP_SIM;

CAchg_SIM              = CAY_SIM(2:end)-CAY_SIM(1:end-1);
CAchg_SIM              = [0;CAchg_SIM];

Lev_SIM                = bpSIM./Y_SIM;  % Leverage
LevSP_SIM              = bpSP_SIM./YSP_SIM;

CtSIM                  = (omega.*cSIM.^(-ita)+(1-omega).*yNSIM.^(-ita)).^(-1/ita);
CtSPSIM                = (omega.*cSP_SIM.^(-ita)+(1-omega).*yNSIM.^(-ita)).^(-1/ita);

RER_SIM                = (omega^(1/(1+ita))+(1-omega)^(1/(1+ita)).*pSIM.^(ita/(1+ita))).^((1+ita)/ita); % Real Exchange Rate
RERSP_SIM              = (omega^(1/(1+ita))+(1-omega)^(1/(1+ita)).*pSP_SIM.^(ita/(1+ita))).^((1+ita)/ita);

