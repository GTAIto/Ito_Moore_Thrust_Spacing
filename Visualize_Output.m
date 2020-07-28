%=========================================================================
%% Plots HIPSTER output
% Called by stepthrough_See_SiStER
%=========================================================================

%------------------------------------------------------------------------
%% Plot PARAMETERS
%% ------------------------------------------------------------------------
fignum=10;
numplots=6;
km=1e3;
MPa=1e6;
vscale=1e3*60*60*24*365;  % for mm/yr
tscale=1e6*60*60*24*365;  % for Myr
Nx = size(X,2); Ny = size(X,1);
dx=diff(x); dy=diff(y);
xn=(x(1:Nx-1)+x(2:Nx))/2;       %coordinates of normal nodes
yn=(y(1:Ny-1)+y(2:Ny))/2;
[XN YN]=meshgrid(xn,yn);
if (exist('Rn','var')==0)
    Rn=zeros(Ny,Nx);
end;
xend=max(x)/1e3;

if exist('seerun','var') ==0
    seerun = ' ';
end;
%% Computes these if numplots>3
if (numplots>3)
    if exist('vxc','var')==0
        SiStER_interp_velocities_to_shear_nodes;
    end;

    %Plot these quivers for velocity
    skipq=6;
    ii=[Ny:-skipq:1]; jj=[1:skipq:Nx];  %starting from the bottom left

    %---------------------------------------------------------------------
    % Stresses (SiStER_output puts them all on normal nodes)

    [n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,sxym);
    SXY  = n2interp(1).data;
    [n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,sxxm);
    SXX  = n2interp(1).data;
    thet=angp(SXX,-SXX,SXY);
    skiptraj=skipq;
    jtr=[1:skiptraj:Nx-1]; itr=[1:skiptraj:Ny-1];

    sigma3=SIG2(SXX(2:Ny,2:Nx),-SXX(2:Ny,2:Nx),SXY(2:Ny,2:Nx))./MPa;
    sigma1=SIG1(SXX(2:Ny,2:Nx),-SXX(2:Ny,2:Nx),SXY(2:Ny,2:Nx))./MPa;
    dsigma=sigma1-sigma3;

    SII=sqrt(SXX.^2 + SXY.^2);

    dY=repmat([0 0.5*dy(1) (dy(2:end)+dy(1:end-1))./2]',1,Nx);
    Plith=cumsum([zeros(1,Nx); [rho(1:Ny-1,1) (rho(1:Ny-1,2:Nx)+rho(1:Ny-1,1:Nx-1))./2];].*dY)*PARAMS.gy;
end
%% Strain rate
if (exist('EXX_s','var'));
    EII=sqrt(0.5*(EXX-Rn./2).^2 + 0.5*(Rn./2-EXX).^2 + EXY.^2);
end

%---------------------------------------------------------------------
%% Failure stress
if (exist('phase_s','var')==0)
    [n2interp]=SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,rhom);
    rhos=n2interp(1).data;
    iwater_s=find(rhos==1000);
    [n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,rhom); 
    rho = n2interp(1).data;
    iwater_n=find(rho(2:end,2:end)==1000);
else
    iwater_s=find(phase_s<2);
    iwater_n=find(phase_n(2:end,2:end)<2);
end
%---------------------------------------------------------------------
if (exist('PARAMS.conv_crit1','var')==0);  %
    PARAMS.conv_crit1=1e-2;
end;

%---------------------------------------------------------------------
%% Stratigraphy: 

[~,imarks]=find(im>1);  %all non-air markers only
strm=stratam(imarks);
if (max(im)==4)
    ibackstop=find(im(imarks)>3);
end
ibase=find(im(imarks) > 2);
ierod=find(strm < 0);
strm=mod(round(strm),2); 
strm=strm+0.28;
strm=min(strm,1);
strm(ierod)=0.75;
strm(ibase)=strm(ibase)-1;

%---------------------------------------------------------------------
%% Theoretical taper angle
if (exist('topo_y','var')==0)
    topo_y=zeros(1,100);
    topo_x=topo_y;
end
phi=atand(MAT(2).mu);
if (length(MAT)<3)
   mubase=PARAMS.decollement_mu;
else
   mubase=MAT(3).mu;
end
bet=atand(abs(PARAMS.gx/PARAMS.gy));

if (isfield(MAT,'lambda_pp')==1)
    lam=MAT(2).lambda_pp; lamb=MAT(3).lambda_pp;
else
    lam=0; lamb=0;
end
rhow=MAT(1).rho0;

if (isfield(PARAMS,'lambda_pp')==1)
    lam=PARAMS.lambda_pp;
    lamb=lam;
else
    lam=0; lamb=0;
end

[alpha,~,psib,~]=critical_wedge_solution(bet,tand(phi),mubase,MAT(2).rho0,MAT(1).rho0,lam,lamb);
taper=alpha+bet;

jtop=min(find(topo_y==min(topo_y))); jbot=max(find(topo_y<(max(topo_y)-100)));
xtaper=[floor(topo_x(jtop)./km):1:ceil(topo_x(jbot)./km)+5]';
if (length(xtaper)==0); xtaper=[1:10]; end
ytaper=min(topo_y)./km+tand(taper).*(xtaper-xtaper(1))-1;


%% ========================================================================
% Plots
%============================================================================
h=0.14; w=0.9; x0=0.08; y0=0.8;ybuf=0.015;
ax=[0 x(end) 0 y(end)]/1e3; markersize=2;
markertype='.';

eval(['figure(' num2str(fignum) '); clf;']) 
set(gcf,'Color','w'); 

%-------------------------------------------------------------------------
% STRATIGRAPHY
%-------------------------------------------------------------------------
subplot('Position',[x0,y0,w,h]);
fastscatter(xm(imarks)/km,ym(imarks)/km,strm,markertype,'markersize',markersize); colorbar; hold on;
caxis([-1 1]);
plot(topo_x/km,topo_y/km,'k--','Linewidth',2)
plot([0 max(topo_x)/km],GEOM(2).top*[1 1]/1e3,'b--');
plot(xtaper,ytaper,'k-')
ylabel('Stratigraphy');     
junk=find(seerun=='_'); modelname=seerun; modelname(junk)=' ';
title([modelname ': time=' num2str(time/tscale) 'Myr; timestep=' num2str(t)],'Fontsize',10);
colormap('jet');
axis ij;
axis equal
axis([0 xend min(y/km) max(y/km)]);
%-------------------------------------------------------------------------
% Remaining plots
%-------------------------------------------------------------------------
vxc(iwater_s)=nan;
thet(iwater_s)=nan;
H=GEOM(2).bot-GEOM(2).top;
pscale=9.8*(MAT(2).rho0*H+MAT(1).rho0.*GEOM(2).top)/1e6;
phi=atand(MAT(2).mu);
yscale=pscale*sind(phi)+MAT(2).Cmax*cosd(phi)/1e6;
for iplot=2:numplots
    clear clim
    if (iplot==2)
        if (exist('epsII_n','var'))
            epsII=epsII_n(2:end,2:end);
        else
             epsII=epsII_s;           
        end
        fieldstr='(epsII)'; clim=[0 0.25].*max(epsII(:));
    elseif (iplot==3)
         fieldstr='(SXX-p)/MPa'; clim=[-2 0]*100; %clim=[-5 0];
    elseif (iplot==4)
        fieldstr='log10(etas)';
	elseif (iplot==5)
        fieldstr='SXY/MPa'; clim=20*[-1 1]; 

    elseif (iplot==6)
        fieldstr='p(2:end,2:end)/MPa'; clim=([-0.1 10]).*pscale/20;
    end

   
    subplot('Position',[x0,y0-(iplot-1).*(ybuf+h),w,h]);
    eval(['field=' fieldstr ';']);
    if (exist('clim','var')==0);
        clim=[min(field(:)) max(field(:))];
    end
    
    if (size(field,1)==Ny-1);
        field(iwater_n)=nan;  
        pcolor(xn/km,yn/km,field); shading flat; hold on;  
     else
        field(iwater_s)=nan;   
        pcolor(x/km,y/km,field); shading flat; hold on;
    end
    grid on
    plot(topo_x/km,topo_y/km,'k-','Linewidth',1)
    if (exist('epsII_n','var'))
        contour(xn./km, yn/km, epsII,0.1*max(epsII(:)),'w');
    else
        contour(x./km, y/km, epsII,0.1*max(epsII(:)),'w');
    end
    
    plot(xtaper,ytaper,'k-')
    ylabel(fieldstr);
    colorbar; 
    caxis(clim); 
    axis ij; 
    axis equal;
    if (iplot==4); quiver(x(jtr)/km,y(itr)/km,vxc(itr,jtr),vyc(itr,jtr),0.75,'k'); end;
    
end
