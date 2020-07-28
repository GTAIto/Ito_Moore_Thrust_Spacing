%=========================================================================
function [alpha1,alpha2,psib1,psib2]= critical_wedge_solution(beta,mu,mub,rho,rhow,lam,lamb)
%Uses Eq.(9), (17), and (19) of Dahlen 1984 to solve for
%alpha, psi0, & psib given mu, mub, lambda, lambdab, rho, & beta
%G.Ito 5/14/17
%=========================================================================
idebug=0;
if (idebug);
    clear;
    rhow=0;
    mu=tand(20); mub=tand(10); rho=2500; lam=0; lamb=0; beta=5;
    mu=1.1; 
    idebug=1;
end;

phi=atand(mu);
phib=atand(mub);

phibp=atand(tand(phib).*(1-lamb)./(1-lam));
mubp=tand(phibp);

psi=[-89:1:89];

fpsi=tand(2.*psi)./(cscd(phi).*secd(2.*psi)-1);
psi(isnan(fpsi))=[];
fpsi(isnan(fpsi))=[];

alph=atand(fpsi*(1-lam)./(1-rhow/rho));
%alphp=alph*(1-rhow./rho)./(1-lam);

psib1=0.5*[asind(sind(phibp)./sind(phi))-phibp];
psib2=0.5*(180-asind(sind(phibp)./sind(phi))-phibp);

bet1=psib1-psi-alph;
bet2=psib2-psi-alph;

alpha1=interp1(bet1,alph,beta);
alpha2=interp1(bet2,alph,beta);

%% ----------------------------------------------------------------
if (idebug==1);
figure(2); clf;
subplot(221);
plot(psi,fpsi); hold on;
plot([psib1 psib2], mubp*[1 1],'o');
xlabel('\psi');
ylabel('F(\psi)')
%xlim([-90 90]);
grid on;

subplot(212); hold on;
plot(bet1, alph,'b'); hold on;
plot(bet2, alph,'b--');
plot([beta beta],[alpha1 alpha2],'o');
ylabel('\alpha');
xlabel('\beta');
%axis([-10 90 -35 35 ]);
%xlim([-10 90]);
grid on;
end;
% ii=find(abs(bet2)==min(abs(bet2)));

