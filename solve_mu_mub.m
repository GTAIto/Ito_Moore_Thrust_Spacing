%=====================================================================
function [phisol, phibsol,phivec,phi_ab,phibw,alphgrid,P1grid,P1calc]=solve_mu_mub(w0,H0,alph_d,bet_d,lam,lamb,coh,rhow,rho,Pparams,nomessage)
%-=====================================================================

rg=rho*9.8;
P1=Pparams(1); P2=Pparams(2); P3=Pparams(3); P4=Pparams(4); yf0=Pparams(5);
dp=0.5;
phivec=[dp:dp:60]'; nn=length(phivec);

phi_ab=NaN*ones(size(phivec));
alphgrid=NaN*ones(nn,nn);
P1grid=alphgrid;
w_d=median(w0)./median(H0);
HH=median(H0)*1e3;
%% ---------------------------------------------------------------
% For each value of mub/phib find the value of phi that gives the observed
% taper (alpha_d+bet_d). i.e., find phi_ab given phib (=phivec). 
% (done by computing grid of alpha for each phib,phi pair and interpolating)
for j = 1:nn  %loop over phib
    mub=tand(phivec(j));
    for i=1:nn  %loop over phi
        mu=tand(phivec(i));
        if (mu>mub)
            [alph1,alph2,psib1,psib2]= critical_wedge_solution(bet_d,mu,mub,rho,rhow,lam,lamb);
            alphgrid(i,j)=alph1;
        end
        FF=(1+sind(phivec(i)))./(1-sind(phivec(i)));
        CCstar=2.*coh./rg.*sqrt(FF);
        wwscale=(FF./2+CCstar./HH+yf0).^P2./((1-lam).*mub+P4*tand(bet_d)).^P3;
        P1grid(i,j)=w_d./wwscale;

    end
    %Finding the values of phi vs phib where alphagrid=alpha_d
    ii=find(~isnan(alphgrid(:,j)));
    if (~isempty(ii))
        if (alph_d>=min(alphgrid(ii,j)) && alph_d <=max(alphgrid(ii,j)))
            phi_ab(j)=interp1(alphgrid(ii,j),phivec(ii),alph_d);
        end
    end
end
%% ------------------------------------------------------------
F=(1+sind(phi_ab))./(1-sind(phi_ab)); Cstar=2.*coh./rg.*sqrt(F);
wsnum=(F./2+Cstar./HH+yf0);
phibw=(((P1.*wsnum.^P2)./w_d).^(1./P3))./(1-lam) - P4*tand(bet_d)./(1-lam);
phibw=atand(phibw);
    
wscale=(F./2+Cstar./HH+yf0).^P2./((1-lam).*tand(phivec)+P4.*tand(bet_d)).^P3;  %scaling for each value of phib, phi that satisfies critical taper
P1calc=w_d./wscale;  %vector of possible values of P1


icr=find(phibw(2:end)>=phivec(2:end) & phibw(1:end-1)<=phivec(1:end-1));
if (length(icr)>0) %Interpolate to find more accurate value
    dphi=phi_ab(icr+1)-phi_ab(icr);
    dpcdp=(phivec(icr+1)-phivec(icr))./dphi;
    dpwdp=(phibw(icr+1)-phibw(icr))./dphi;
    deltaphi=(phibw(icr)-phivec(icr))./(dpcdp-dpwdp);
    phibsol=phivec(icr)+deltaphi*dpcdp;
    phisol=phi_ab(icr)+deltaphi;
    nosol=0;
   
else
    icr=find(abs(phibw-phivec)==min(abs(phibw-phivec)),1);
    phibsol=(phibw(icr)+phivec(icr))/2;
    phisol=phi_ab(icr);
    num1=sqrt((phibw(icr)-phivec(icr)).^2);
    den1=sqrt(phibsol.^2+phi_ab(icr).^2);
    normerror=num1./den1;
    nosol=1;
    P1_error=min(abs(P1-P1calc))./P1;
    if (nomessage==0); disp(['>>No exact solution, closest values are output. Error=' num2str(num1), ' Normerror=' num2str(normerror) ' P1_error=' num2str(P1_error)]); end;
end

F=(1+sind(phi_ab))./(1-sind(phi_ab)); Cstar=2.*coh./rg.*sqrt(F);
wscale=(F./2+Cstar./HH+yf0).^P2./((1-lam).*tand(phivec)+P4.*tand(bet_d)).^P3;  %scaling for each value of phib, phi that satisfies critical taper
P1test=w_d./wscale;


if (nosol==1);
    ii=find(P1test<=P1,1);
    if (isempty(ii))
        ii=find(P1test==min(P1test));
    end
        
    phibsol=phivec(ii);
    phisol=phi_ab(ii);
end













