%=============================================================================
% Updates stresses on markers for CURRENT solution.  Stress rotation
% occurs after solutions are output
% G.Ito 8/16
%==============================================================================

% Compute STRESS Changes, interpolate to markers, and apply to marker stresses
% Remember: SYY=-SXX because these are stress deviators
dsxx=(2*etan.*EXX-etan.*Rn - sxxOLD).*Zn;
dsxy=(2*etas.*EXY-sxyOLD).*Zs;

[dum]=SiStER_interp_normal_nodes_to_markers(dsxx,xc,yc,xm,ym,icn,jcn);
sxxm=sxxm+dum;

[dum]=SiStER_interp_shear_nodes_to_markers(dsxy,x,y,xm,ym,icn,jcn);
sxym=sxym+dum;

%--------------------------------------------------------------------------------
% Pressures
% Added G.Ito 5/18
%--------------------------------------------------------------------------------
[dum]=SiStER_interp_normal_nodes_to_markers((p-pold),xc,yc,xm,ym,icn,jcn);
pm=pm+dum;