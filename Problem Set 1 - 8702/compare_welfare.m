%% 5. Welfare Calculation
V    = zeros(NB,NSS);
VSP  = zeros(NB,NSS);

U    = (totalcDE.^(-1/ita)).^(1-sigma)./(1-sigma);
USP  = (totalcSP.^(-1/ita)).^(1-sigma)./(1-sigma);

EV   = U;
EVSP = USP;

d3   = 100;
iter = 0;
disp('Value Iter     Norm');
while d3>tol && iter <iter_tol
    for i=1:NB
        for j=1:NSS
            EV(i,j)   = beta*interp1(B,V,bpDE(i,j),'linear','extrap')*Prob(j,:)'; %EV is expected value tomorrow in today's grid.
            EVSP(i,j) = beta*interp1(B,VSP,bpSP(i,j),'linear','extrap')*Prob(j,:)';
        end
    end
    
    V_new   = U+EV;
    VSP_new = USP+EVSP;
    
    d3      = max([max(max(abs((V_new-V)./V))),max(max(abs((VSP_new-VSP)./VSP)))]);
    
    iter    = iter+1;
    
    V       = V_new;
    VSP     = VSP_new;
    if mod(iter, outfreq) == 0
        fprintf('%d          %1.7f \n',iter,d3);
    end
    
end
Wf=(VSP./V).^(1/(1-sigma))-1;

%% 6. Optimal Tax
 
mupSP    = (omega*totalcSP.^(sigma/ita).*(totalcSP.^(-1/ita-1)).*(cSP.^(-ita-1))); %mupSP-EESP.*psi;
for i=1:NB
    for j=1:NSS
        emuSP(i,j) = beta*(SR(i,j))*interp1(B,mupSP,bpSP(i,j),'linear','extrap')*Prob(j,:)';
    end
end
tau           = mupSP./emuSP-1;


% tau2= 

tau(IbindSP ==1) = 0;
tau(tau<0)    = 0;

% figure('name','tax as a function of debt')
% plot(B,tau); legend('1','2','3','4','5','6','7','8','9')
% xlabel('debt')
% figure('name','tax as a function of income')
% for IB=1:10:50
%     plot(yT,tau(IB,:)); xlabel('income'); title('tax');     hold on
% end
% 
% plots_comp
