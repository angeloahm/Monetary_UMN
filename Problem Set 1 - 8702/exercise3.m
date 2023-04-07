%% 1-2. Model/numerical Parameter Values
% PARAMETERS
sep_prefs=1;

sigma    = 2;
ita      = (1/0.83)-1; % elasticity parameter (elasticity is 0.83)
beta     = 0.906;       % discount factor
omega    = 0.3070;       % share for tradables
kappa    = 0.3235;

R_grid = [1.03 1.05];
yT = [0.91 1.02];
yN = yT;
Pi_R = [0.8 0.2; 0.2 0.8];
Pi_y = [0.2 0.1; 0.8 0.9]';

bmin_NDL    = -kappa*(1+min(yT(1)))+0.01;  % maximum debt consistent with feasibility
bmin        = -1.02;
bmax        = -0.2000;


NB          = 100;
B           = linspace(bmin,bmax,NB)';


SR          = ones(NB,2,2);
SR(:, :, 1) = R_grid(1);
SR(:, :, 2) = R_grid(2);

YT          = ones(NB,2,2);
YT(:, 1,:)  = yT(1);
YT(:, 2,:)  = yT(2);

YN          = ones(NB,2,2);
YN(:, 1,:)  = yN(1);
YN(:, 2,:)  = yN(2);

b           = ones(NB,2,2);
for i=1:NB
    b(i, :, :) = B(i);
end

KAPPAS      = kappa*ones(NB,2,2);

% ininitialize
bp          = b;
c           = b  +YT  -bp./SR;
 


price       = (1-omega)/omega*c.^(1+ita);
EE          = zeros(NB,2,2);
bmax_collat = -kappa*(price.*YN +YT);
cbind       = b  +YT  -bmax_collat./SR;
emu         = zeros(NB,2,2)*nan;
Ibind       = ones(NB,2,2);
 
 
uptd        = 0.2;   % Change to 0.1 if policy function fluctuate
uptdSP      = 0.02;

outfreq     = 50;    % Display frequency
iter_tol    = 500;   % Iteration tolerance
Tol_eulerDE = 1.0e-7;% when to look for unconstrained solution
Tol_eulerSP = 1.0e-7;
tol         = 1.0e-5; %Numerical tolerance


cut = 0;  % Initial condition dependence


options = optimset('Display','off');

disp(' ')


%% 3) Decentralized Equilibrium

iter = 0;
d2   = 100;

tic

%figure('name','bp')

disp('DE Iter      Norm         Updt');

while d2>tol && iter <iter_tol
    
    oldp   = price;
    oldc   = c;
    oldbp  = bp;
    
    totalc = (omega*c.^(-ita)+(1-omega)*YN.^(-ita));
    mup    = real(omega*totalc.^(sigma/ita).*(totalc.^(-1/ita-1)).*(c.^(-ita-1)));  % marginal utility today
    
          if sep_prefs==1
                 mup    = omega*c.^(-sigma);
          end
    
    for i=1:NB
        for j=1:2
            for k=1:2
%                 if i>1
%                     keyboard
%                 end
                mup_int=interp1(B,mup,bp(i,j,k),'linear','extrap');
                ER_mup = sum(mup_int(1,j,:) .* Pi_R(k,:));
                Ey_mup = sum(mup_int(1,:,k) .* Pi_y(j,:));
                emu(i,j,k) = beta*sum(Ey_mup(:)+ER_mup(:)); %EMU is expected marginal utility tomorrow in today's grid.
            end
        end
    end
    
    
    for i=1:NB
        for j=1:2
            for k=1:2
                EE(i,j,k) = (omega*cbind(i,j,k)^(-ita)+(1-omega)*YN(i,j,k)^(-ita))^(sigma/ita -1/ita -1 )*omega*cbind(i,j,k)^(-ita-1)-emu(i,j,k); 
            end
            
            if EE(i,j,k)>Tol_eulerDE     % BORROWING CONSTRAINT BINDS
                bp(i,j,k)    = bmax_collat(i,j,k);
                c(i,j,k)     = cbind(i,j,k);
                Ibind(i,j,k) = 1; 
      
            else                        % EULER EQUATION
                
                
                if sep_prefs ==0   
                    f                = @(cc) (omega*cc^(-ita)+(1-omega)*YN(i,j,k)^(-ita))^(sigma/ita -1/ita -1)*omega*cc^(-ita-1)-emu(i,j,k);
                    c0               = c(i,j,k);
 
                    [c(i,j,k),EE(i,j,k)] = fzero(f,c0,options);
                else
                    c(i,j,k)           = (emu(i,j,k)/omega) ^(-1/sigma);
                    
                    EE(i,j,k)          = 0;
                end
               
                Ibind(i,j,k)       = 0;
            end
                  
            
        end
    end
    
    bp    = SR.*(b+YT-c);
    bp    = min(bmax,max(bp,bmin));
    
    price = (1-omega)/omega*(c./YN).^(1+ita);
    
    c     = SR.*(b+YT-max(bp,-KAPPAS.*(price.*YN+YT))); % update consumption based on update for pN
    
    d2    = max([max(max(abs(c(:)-oldc(:)))),max(max(abs(bp(:)-oldbp(:))))]); % metric
    
    bmax_collat                   = -KAPPAS.*SR.*(price.*YN+YT);
    bmax_collat(bmax_collat>bmax) = bmax;
    bmax_collat(bmax_collat<bmin) = bmin;
    cbind                         = SR.*(b+YT - bmax_collat );
    cbind(cbind<0)                = nan;
    
    
    %=====================Updating rule.Must be slow, important==============
    bp    = uptd*bp+(1-uptd)*oldbp;
    c     = uptd*c+(1-uptd)*oldc;
 
    %========================================================================
    

    iter              = iter+1;
    
    D2(iter)          = d2;
    
    if mod(iter, outfreq)==0 | iter==1
     fprintf('%3i          %10.2e   %2.2f \n',iter,d2, uptd);
     %plot(B,bp); xlabel('B'); ylabel('B''');  drawnow
    end
end

     fprintf('No. Iterations: %3i  Metric  %10.2e \n',iter,d2 );

LagrangeDE            = EE;
LagrangeDE(Ibind ==0) = 0;

%plots_de
 
bpDE             = bp;
cDE              = c;
totalcDE         = totalc;
priceDE          = price;
IbindDE          = Ibind;
clear bp c totalc  price mup Ibind
%% Plots for the DE
figure(1)
plot(B, cDE(:,1,1), B, cDE(:,1,2), B, cDE(:,2,1), B, cDE(:,2,2), LineWidth=2)
legend('low y - low R', 'low y - high R', 'high y - low R', 'high y - high R')
title('Consumption')
saveas(gcf, 'consumption3.png')

figure(2)
plot(B, bpDE(:,1,1), B, bpDE(:,1,2), B, bpDE(:,2,1), B, bpDE(:,2,2), LineWidth=2)
legend('low y - low R', 'low y - high R', 'high y - low R', 'high y - high R')
title('Debt')
saveas(gcf, 'debt3.png')

%% 4) Planner problem

bp         = bpDE;
c          = cDE;
psi        = (1-omega)/omega*(ita+1)*KAPPAS.*(c./YN).^(ita);
price      = priceDE;
LagrangeSP = LagrangeDE./(1-psi);  % guess for SP
IbindSP    = IbindDE;

EESP       = zeros(NB,2,2)*nan;

uptd = uptdSP;
d2   = 100;
iter = 0;
disp('SPP Iter      Norm       Uptd');

while d2>tol && iter <iter_tol*2
    
 
    oldc          = c;
    oldbp         = bp;
    oldLagrangeSP = LagrangeSP;
    
    totalc        = (omega*c.^(-ita)+(1-omega)*YN.^(-ita));
    mup           = (omega*totalc.^(sigma/ita).*(totalc.^(-1/ita-1)).*(c.^(-ita-1)));  %mup(b,y)
    if ~isreal(mup(:))
        keyboard
    end
          if sep_prefs==1
                 mup    = omega*c.^(-sigma);
          end
    
    psi      = (1-omega)/omega*(ita+1)*KAPPAS.*(c./YN).^(ita);
    LFH_SP   = mup+LagrangeSP.*psi;
    
    for i=1:NB
        for j=1:2
            for k=1:2
                int_LHS=interp1(B,LFH_SP,bp(i,j,k),'linear','extrap');
                ER_LHS = sum(int_LHS(1,j,:) .* Pi_R(k,:));
                Ey_LHS = sum(int_LHS(1,:,k) .* Pi_y(j,:));
                RHS_SP(i,j,k) = beta*sum(ER_LHS(:)+Ey_LHS(:)); %Expected marginal utility tomorrow in today's grid.
            end
        end
    end
    
    LagrangeSP = LFH_SP-beta*RHS_SP; % THIS IS EQUATION u_T +mu*psi= beta*R E [(u_T(t+1)+mu_t+1 ]+mu_t
    
    LagrangeSP(IbindSP ==0) = 0;
    
    for i=1:NB
        for j=1:2
            for k=1:2
                if LagrangeSP(i,j,k)>=Tol_eulerSP          % BORROWING CONSTRAINT BINDS
        
                    c(i,j,k) = cbind(i,j,k);
                    IbindSP (i,j,k) = 1;
                    EESP(i,j,k)= LagrangeSP(i,j,k);
                else
                 
                    if sep_prefs==0                        % EULER EQUATION   
                        f          = @(cc) (omega*cc^(-ita)+(1-omega)*YN(i,j)^(-ita))^(sigma/ita-1/ita-1)*omega*cc^(-ita-1)-RHS_SP(i,j,k);
                        x0         = c(i,j,k);
                        [c(i,j,k),EESP(i,j,k)] = fzero(f,x0,options);
                    else
                        c(i,j,k)=   (RHS_SP(i,j,k)/omega) ^(-1/sigma);
                        EESP(i,j,k)  = 0;
                    end
  
                    IbindSP (i,j,k) = 0;
                end
            end
        end
    end
    
    bp          = SR.*(b+YT-c);
    bp(bp>bmax) = bmax;
    bp(bp<bmin) = bmin;
    price       = ((1-omega)/omega*(c./YN).^(1+ita));
    
    %===============Check collateral constraint==============================
    c           = SR.*(b+YT-max(bp,-KAPPAS.*(price.*YN+YT)));
    price       = ((1-omega)/omega*(c./YN).^(1+ita));
    bp          = SR.*(b+YT-c);
    %========================================================================
 
    d2          = max(max(abs(c(:)-oldc(:))));
  
 
    
    %=====================Updating rule.Must be slow, important==============
    bp          = uptd*bp+(1-uptd)*oldbp;
    c           = uptd*c+(1-uptd)*oldc;
    %========================================================================
    
    bmax_collat = -KAPPAS.*(price.*YN+YT);
    cbind       = SR.*(b+YT-bmax_collat);
    
    iter=iter+1;
   if mod(iter, outfreq) == 0 | iter==1
        fprintf('%3i          %10.2e   %2.2f \n',iter,d2, uptd);
      
   end
 
end

    fprintf('No. Iterations: %3i  Metric  %10.2e \n ',iter,d2);


toc
totalcSP = totalc;
 mupSP    = mup;
 bpSP     = bp;
 priceSP  = price;
 cSP      = c;

clear bp c totalc  price mup

