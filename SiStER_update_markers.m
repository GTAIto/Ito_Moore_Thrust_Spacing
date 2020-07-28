% SiStER Update Markers
% Added stratam for passively visualizing deformation G.Ito 7/18
% Added pm G.Ito 5/18


[xm_new,ym_new] = SiStER_advect_markers(x,y,xm,ym,dx,dy,dt_m,vx,vy);
xm=xm_new;
ym=ym_new;
% eliminate markers that left domain
Iin=find(xm<=xsize & xm>=0 & ym>=0 & ym<=ysize);

msg2='  markers removed: ';
msg=[msg2 num2str(length(xm)-length(Iin))];
disp(msg)
xm=xm(Iin);
ym=ym(Iin);
im=im(Iin);
ep=ep(Iin);
Tm=Tm(Iin);
idm=idm(Iin);
sxxm=sxxm(Iin);
sxym=sxym(Iin);
epsIIm=epsIIm(Iin);
stratam=stratam(Iin);
pm=pm(Iin);  %5/18


% locate advected markers with respect to the eulerian grid
[quad,icn,jcn] = SiStER_locate_markers_in_grid(xm,ym,x,y,dx,dy);
    
% check for holes in the marker distribution, 
% patch with new markers if necessary
% those new markers immediately get assigned a value of phase (im), index 
% (idm) and accumulated plastic strain (ep), i.e., the 2 variables that never get
% passed to nodes. 
[xm, ym, im, Ifix, mp, ep, idm, Tm, sxxm, sxym, pm, epsIIm, stratam]=...
    SiStER_patch_marker_holes(icn,jcn,quad,Nx,Ny,Mquad,Mquad_crit,xm,ym,x,y,dx,dy,im,ep,idm,Tm,sxxm,sxym,pm,epsIIm,stratam);   %G.Ito

% then they get assigned P, epsII and stresses from grid values

if min(Ifix)>0
    
    msg2='  markers added: ';
    msg=[msg2 num2str(length(Ifix))];
    disp(msg)
    xmFIX=xm(Ifix);
    ymFIX=ym(Ifix);
    
    % pass temperature, pressure, strain rate and stresses to the new
    % markers from their nodal values
    
    % locate new markers with respect to the eulerian grid
    [quadFIX,icnFIX,jcnFIX] = SiStER_locate_markers_in_grid(xmFIX,ymFIX,x,y,dx,dy);
    
    
    % interpolate grid values to the new markers
    % from shear nodes:

    
     % dont do that- you'd be passing pre-advection nodal Ts to markers that were just advected    
     %[temp]=SiStER_interp_shear_nodes_to_markers(T,x,y,xmFIX,ymFIX,icnFIX,jcnFIX);
     %Tm(Ifix)=temp; % temperature
     
% UNCOMMENT WHEN I ADD  STRESSES    
%     [temp]=SiStER_interp_shear_nodes_to_markers(sxy0,x,y,xmFIX,ymFIX,icnFIX,jcnFIX);
%     sxym(Ifix)=temp; % stress_xy
%     
%     % from normal nodes
% 
%     [temp]=SiStER_interp_normal_nodes_to_markers(sxx0,xc,yc,xmFIX,ymFIX,icnFIX,jcnFIX);
%     sxxm(Ifix)=temp; % stress_xx
% removed G.Ito 5/4/18
%     [temp]=SiStER_interp_normal_nodes_to_markers(p,xc,yc,xmFIX,ymFIX,icnFIX,jcnFIX);
%     pm(Ifix)=temp; % pressure
    
    % DO IT ALSO FOR STRAIN RATE !
    
end
    
    
% locate all markers with respect to the eulerian grid

[qd,icn,jcn] = SiStER_locate_markers_in_grid(xm,ym,x,y,dx,dy);