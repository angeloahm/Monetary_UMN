%% 4) Planner problem

bp         = bpDE;
c          = cDE;
psi        = (1-omega)/omega*(ita+1)*KAPPAS.*(c./YN).^(ita);
price      = priceDE;
LagrangeSP = LagrangeDE./(1-psi);  % guess for SP
IbindSP    = IbindDE;

EESP       = zeros(NB,NSS)*nan;

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
    
          if sep_prefs==1
                 mup    = omega*c.^(-sigma);
          end
    
    psi      = (1-omega)/omega*(ita+1)*KAPPAS.*(c./YN).^(ita);
    LFH_SP   = mup+LagrangeSP.*psi;
    
    for i=1:NB
        for j=1:NSS
            RHS_SP(i,j) = beta*SR(i,j)*interp1(B,LFH_SP,bp(i,j),'linear','extrap')*Prob(j,:)'; %Expected marginal utility tomorrow in today's grid.
        end
    end
    
    LagrangeSP = LFH_SP-beta*SR.*RHS_SP; % THIS IS EQUATION u_T +mu*psi= beta*R E [(u_T(t+1)+mu_t+1 ]+mu_t
    
    LagrangeSP(IbindSP ==0) = 0;
    
    for i=1:NB
        for j=1:NSS
       
            if LagrangeSP(i,j)>=Tol_eulerSP          % BORROWING CONSTRAINT BINDS
        
                c(i,j) = cbind(i,j);
                IbindSP (i,j) = 1;
                EESP(i,j)= LagrangeSP(i,j);
            else
                 
                if sep_prefs==0                           % EULER EQUATION   
                    f          = @(cc) (omega*cc^(-ita)+(1-omega)*YN(i,j)^(-ita))^(sigma/ita-1/ita-1)*omega*cc^(-ita-1)-RHS_SP(i,j);
                    x0         = c(i,j);
                    [c(i,j),EESP(i,j)] = fzero(f,c0,options);
                else
                    c(i,j)=   (RHS_SP(i,j)/omega) ^(-1/sigma);
                    EESP(i,j)  = 0;
                end
  
                IbindSP (i,j) = 0;
           
            end
        end
    end
    bp          = (SR.*b+YT-c);
    bp(bp>bmax) = bmax;
    bp(bp<bmin) = bmin;
    price       = ((1-omega)/omega*(c./YN).^(1+ita));
    
    %===============Check collateral constraint==============================
    c           = SR.*b+YT-max(bp,-KAPPAS.*(price.*YN+YT));
    price       = ((1-omega)/omega*(c./YN).^(1+ita));
    bp          = (SR.*b+YT-c);
    %========================================================================
 
    d2          = max(max(abs(c-oldc)));
    
  
 
    
    %=====================Updating rule.Must be slow, important==============
    bp          = uptd*bp+(1-uptd)*oldbp;
    c           = uptd*c+(1-uptd)*oldc;
    %========================================================================
    
    bmax_collat = -KAPPAS.*(price.*YN+YT);
    cbind       = SR.*b+YT-bmax_collat;
    
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

