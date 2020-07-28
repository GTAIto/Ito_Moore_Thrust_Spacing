%=========================================================================
% Produces Fig. 14 of Ito & Moore (2020) for solutions to
% Coefficients of friction for Makran, Nankai, Cascadia, and Hikurangi
% Based on scaling laws and measurements of w0/H and beta 
% AND critical Coulomb wedge theory with measurements of alpha & beta
%-========================================================================
%% ---------------------------------------------------------------

clear;
figure(1); clf;
set(gcf,'Color','w'); 
%% Parameters
rhow=1000; rho=2300; 
lam=0.75;
lamb=lam;
coh=0e6; 
rg=rho*9.8;
P1=1.1; P2=0.64; P3=0.16; P4=4.5; yf0=1/2; 
Pparams=[P1 P2 P3 P4 yf0];

%% Observations
dat=load('Obs_widths.txt');
[m,n]=size(dat);
alph=median(dat(:,1:2),2); 
bet=median(dat(:,3:4),2); 
H1=dat(:,5); H2=dat(:,6);
H=median([H1 H2],2);
w1=dat(:,7); w2=dat(:,8);
w=median([w1 w2],2);
ndso=dat(:,9);

w=w.*cosd(bet);
H=H.*cosd(bet);
Hso=zeros(max(ndso),1);
betso=Hso;
for n=1:max(ndso)
    Hso(n)=median(H(ndso==n))*1e3;
    betso(n)=median(bet(ndso==n));
end
taper=alph+bet;


%% ---------------------------------------------------------------
% Plotting variables
fsize=10;
t=1; th=1;
lw=[t t t t t t t t t th th th th th];

bb=0.2; cc=0.8;

s=11; b=16; 
sz=[s s s s s s s s s b b b b b];
mcol=['ks'; 'ko'; 'k^'; 'kd'; 'ko'; 'ko';];
col=[0.7*[1 1 1]; [bb cc bb]; [bb bb cc]; [cc bb bb]; [cc cc bb]; [cc cc bb];];

n=length(w);

%------------------------------------------------------------
%% PLOTS
h=0.25; ww=0.35;
xs=[0.11 0.49 0.11 0.49 0.49 0.3]; ys=[0.4 0.4 0.12 0.12 0.12 0.7];

colordef white
Hvec=[0:20];
ndsop=[1 2 3 4];
subplot('Position',[xs(end) ys(end) ww h]); 
    for j=1:4
        plot(Hvec,Hvec*j,'k','LineWidth',0.5);  hold on;
        if (j<4)
           plot(Hvec,Hvec*(j+0.5),'k--','LineWidth',0.5);  hold on;
        end
    end

for n=1:length(ndsop)
ndplot=find(ndso==ndsop(n));
kmin=min(ndplot); kmax=max(ndplot); 
subplot('Position',[xs(end) ys(end) ww h]); 
for k=kmin:kmax
    errorbar(H(k),w(k),w(k)-w1(k),w2(k)-w(k),mcol(ndsop(n),:),'MarkerFaceColor',col(ndsop(n),:),'LineWidth',lw(ndsop(n)),'MarkerSize',sz(ndsop(n))); hold on;
    set(gca,'FontSize',fsize)
    axis([0 5 0 14]);
end

%% ------------------------------------------------------------
% Solutions for phi and phib
alphcont=[0:1:7];
w0cont=[0:0.25:4];

for k=kmin:kmax
    bet_d=bet(k);
    alph_d=alph(k);
    w0=w(k);
    H0=H(k);
    subplot('Position',[xs(ndsop(n)) ys(ndsop(n)) ww h]);
    if(k==kmin)
        fprintf('****ndso= %3.0f, Hmedian=%3.2f, betamed=%2.1f \n', [ndso(k) Hso(ndso(k))/1e3 betso(ndso(k))]);
        nomessage=1;
        [phisol,phibsol,phivec,phi_ab,phibw,alphgrid,P1grid,P1calc]=solve_mu_mub(w0,Hso(ndso(k)),alph_d,betso(ndso(k)),lam,lamb,coh,rhow,rho,Pparams,nomessage);
        calph=contour(tand(phivec),tand(phivec),alphgrid,alphcont,'k','LineWidth',0.3); hold on;
        [phibg phig]=meshgrid(phivec,phivec); 
        F=(1+sind(phig))./(1-sind(phig)); Cstar=2.*coh./rg.*sqrt(F);
        w0norm=P1.*((F./2+Cstar./Hso(ndso(k))+yf0).^P2)./((1-lam).*tand(phibg)+P4.*tand(betso(ndso(k)))).^P3;
        cw0=contour(tand(phivec),tand(phivec),w0norm,w0cont,'Color',[0.7 0 0],'Linewidth',0.3); 
        clabel(cw0);
    end
    [phisol1,phibsol1,phivec,phi_ab,phibw,alphgrid,P1grid,P1calc]=solve_mu_mub(w1(k),H0,alph_d,bet_d,lam,lamb,coh,rhow,rho,Pparams,1);
 	[phisol2,phibsol2,phivec,phi_ab,phibw,alphgrid,P1grid,P1calc]=solve_mu_mub(w2(k),H0,alph_d,bet_d,lam,lamb,coh,rhow,rho,Pparams,1);
     nomessage=0;
    [phisol,phibsol,phivec,phi_ab,phibw,alphgrid,P1grid]=solve_mu_mub(w0,H0,alph_d,bet_d,lam,lamb,coh,rhow,rho,Pparams,nomessage);
    %phisol1=0.9*phisol; phisol2=1.1*phisol; phibsol1=0.9*phibsol; phibsol2=1.1*phibsol;
    fprintf('nd=%3.0f alph=%2.1f, bet=%2.1f w0=%3.2f H=%3.2f w0/H=%3.2f mub%4.2f, mu=%4.2f, Hso=%3.2f, betso=%2.1f, ndso=%3.0f\n',....
             [k alph_d bet_d w0 H0 w0./H0 tand(phibsol) tand(phisol) Hso(ndso(k))/1e3 betso(ndso(k)) ndso(k)])

    nnot=find(ndsop~=ndsop(n));
    for nn=nnot
        subplot('Position',[xs(ndsop(nn)) ys(ndsop(nn)) ww h]);
        %errorbar(phibsol,phisol,phisol-phisol1,phisol2-phisol,phibsol-phibsol1,phibsol2-phibsol,mcol(ndsop(n),:),'MarkerFaceColor',[1 1 1],'LineWidth',lw(ndsop(n)),'MarkerSize',sz(ndsop(n))); hold on;
        plot(tand(phibsol),tand(phisol),mcol(ndsop(n),:),'MarkerFaceColor',[1 1 1],'LineWidth',lw(ndsop(n)),'MarkerSize',sz(ndsop(n))); hold on;

    end
    subplot('Position',[xs(ndsop(n)) ys(ndsop(n)) ww h]);
    errorbar(tand(phibsol),tand(phisol),tand(phisol-phisol1),tand(phisol2-phisol),tand(phibsol-phibsol1),tand(phibsol2-phibsol),....
        mcol(ndsop(n),:),'MarkerFaceColor',col(ndsop(n),:),'LineWidth',lw(ndsop(n)),'MarkerSize',sz(ndsop(n))); hold on;
    %Draw boxes for lab-derived solutions
    if (ndso(k)==2)  
        mb1=0.27; mb2=0.31; mw1=0.30; mw2=0.38;
        plot([mb1 mb1 mb2 mb2 mb1],[mw1 mw2 mw2 mw1 mw1],'k'); hold on;
        mb1=0.28; mb2=0.31; mw1=0.24; mw2=0.38;
        plot([mb1 mb1 mb2 mb2 mb1],[mw1 mw2 mw2 mw1 mw1],'k')
        mb1=0.2; mb2=0.25; mw1=0.23; mw2=0.28;
        plot([mb1 mb1 mb2 mb2 mb1],[mw1 mw2 mw2 mw1 mw1],'g')
    elseif (ndso(k)==3)
        mb1=0.18; mb2=0.20; mw1=0.19; mw2=0.75;
        plot([mb1 mb1 mb2 mb2 mb1],[mw1 mw2 mw2 mw1 mw1],'b')
    elseif (ndso(k)==4);
        mb1=0.23; mb2=0.39; mw1=0.33; mw2=0.48;
        mb1=0.07; mb2=0.22; mw1=0.2; mw2=0.35;
        plot([mb1 mb1 mb2 mb2 mb1],[mw1 mw2 mw2 mw1 mw1],'r')

    end

end
axis([0 0.8 0 1.5])
set(gca,'FontSize',fsize)

end

%export_fig -pdf -r600 Fig14_Obs_matlab