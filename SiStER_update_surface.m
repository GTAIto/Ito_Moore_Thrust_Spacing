%% ========================================================================
% SiStER_update_surface
% advects surface, erodes it, and then ensures they remain across the
% surface for contractional cases
% G. Ito modified from J. Weiss's code 11/15
% G.Ito overhanging topo. now removed BEFORE erosion calculation 7/4/17
%% ========================================================================

if (exist('isurface','var'))
if (isurface==1)
    if (BC.right(3)>=0);
        disp('Surface topography works only for inflow on left side');
        disp('Stopped in SiStER_update_surface.m')
        halt
    end
    % UPDATE MARKER CHAIN FOR TOPOGRAPHY
    [topo_x,topo_y] = SiStER_advect_markers(x,y,topo_x,topo_y,dx,dy,dt_m,vx,vy);
    %Remove any cliffs or overhanging topography
    ii=find(topo_x(2:end)-topo_x(1:end-1)<(xsize./Ntopo/1e3));
    if (length(ii)>=1);
       topo_x(ii+1)=[];
       topo_y(ii+1)=[];
    end;
    %topo_x must remain in order of increasing x. 
    [topo_x,isort]=sort(topo_x);
    topo_y=topo_y(isort);
% MAKE SURE TOPO MARKERS SPAN FULL WIDTH OF DOMAIN
% Original code for extension case
%    Iin=find(topo_x<=xsize & topo_x>=0);
%    topo_x=topo_x(min(Iin)-1:max(Iin)+1);
%    topo_y=topo_y(min(Iin)-1:max(Iin)+1);
    num=round((xsize-topo_x(end))/PARAMS.dxsurf); % number of markers to add (~10 m increment)
    if (num>0)
        new_x=linspace(topo_x(end)+PARAMS.dxsurf,xsize,num); % new x values
        new_y=topo_y(end)*ones(size(new_x)); % new y values
        topo_x=[topo_x new_x];
        topo_y=[topo_y new_y];
    end;

% Reseed topography if needed
    reseed_yes=0;
    if (PARAMS.kappa>0);
        dt_surf=0.5*min(diff(topo_x)).^2/PARAMS.kappa;
        nsolve=ceil(dt_m/dt_surf);
        if (nsolve>1e3);
            reseed_yes=1;
        end;
    end;
    if (mod(t,PARAMS.reseedsurface)==0 || reseed_yes)
        disp(['Reseeding topography t = ' num2str(t)])
        sorter = sortrows([topo_x' topo_y'],1);    % makes sure things are in the ride order
        topo_x = sorter(:,1)'; 
        topo_y = sorter(:,2)';
        %Ntopo  = 1e5;   % number of topo markers
        Ntopo=round(xsize/PARAMS.dxsurf)+1;
        dum1  = topo_x;
        topo_x =linspace(0,xsize,Ntopo);
        topo_y = interp1(dum1,topo_y,topo_x);
    end
% EROSION
tcpu4=cputime;
    if PARAMS.kappa>0    
        topo_y=erosion_by_diffusion(topo_x,topo_y,dt_surf,dt_m,PARAMS.kappa);
    
    % % RESET TOPOGRAPHY
    % % i.e. make everything above/below surface air/rock
    % % added by GI/JW/JAO - 07/15
        %ind=find(isnan(topo_y));topo_x(ind)=[];topo_y(ind)=[]; % interp1 choking on NaN/Inf values for some diffusivities - JW 07/15
        ii=find(topo_x>-0.1*xsize & topo_x<1.1*xsize & topo_y>-0.1*ysize & topo_y<1.1*ysize);
        topo_x=topo_x(ii); topo_y=topo_y(ii);
      %Removing overhanging topomarkers WAS here
        topomarkers=interp1(topo_x,topo_y,xm,'spline','extrap');
        im(ym<topomarkers & im==2)=1;
        ii=ym>=topomarkers & im==1;
        im(ii)=2;
        stratam(ii)=-2;
        
    end
end
end
tcpu5=cputime;
