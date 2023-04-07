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
        for j=1:NSS
            emu(i,j) = beta*(SR(i,j)*STAU(i,j))*interp1(B,mup,bp(i,j),'linear','extrap')*Prob(j,:)'; %EMU is expected marginal utility tomorrow in today's grid.
 
        end
    end
    
    
    for i=1:NB
        for j=1:NSS
            EE(i,j) = (omega*cbind(i,j)^(-ita)+(1-omega)*YN(i,j)^(-ita))^(sigma/ita -1/ita -1 )*omega*cbind(i,j)^(-ita-1)-emu(i,j); 
           
            
            if EE(i,j)>Tol_eulerDE     % BORROWING CONSTRAINT BINDS
                bp(i,j)    = bmax_collat(i,j);
                c(i,j)     = cbind(i,j);
                Ibind(i,j) = 1; 
      
            else                        % EULER EQUATION
                
                
                if sep_prefs ==0   
                    f                = @(cc) (omega*cc^(-ita)+(1-omega)*YN(i,j)^(-ita))^(sigma/ita -1/ita -1)*omega*cc^(-ita-1)-emu(i,j);
                    c0               = c(i,j);
 
                    [c(i,j),EE(i,j)] = fzero(f,c0,options);
                else
                    c(i,j)           = (emu(i,j)/omega) ^(-1/sigma);
                    
                    EE(i,j)          = 0;
                end
               
                Ibind(i,j)       = 0;
            end
                  
            
        end
    end
    
    % (ADD THE REBATE HERE!)
    bp    = STAU.*SR.*b+YT-c+(1-STAU).*b;
    bp    = min(bmax,max(bp,bmin));
    
    price = (1-omega)/omega*(c./YN).^(1+ita);
    
    c     = STAU.*SR.*b+YT-max(bp,-KAPPAS.*(price.*YN+YT))+(1-STAU).*b; % update consumption based on update for pN
    % (ADD TAU IN THIS EQUATION!)
    
    d2    = max([max(max(abs(c-oldc))),max(max(abs(bp-oldbp)))]); % metric
    
    bmax_collat                   = -KAPPAS.*(price.*YN+YT);
    bmax_collat(bmax_collat>bmax) = bmax;
    bmax_collat(bmax_collat<bmin) = bmin;
    cbind                         = STAU.*SR.*b+YT - bmax_collat + (1-STAU).*b;
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