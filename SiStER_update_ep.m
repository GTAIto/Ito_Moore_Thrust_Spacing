
%=========================================================================
% PLASTIC STRAIN ACCUMULATION
% 
% G.Ito 8/16  original version using epsII_s from flow slove, which involves 
%             interpolating EXX to shear nodes
% G.Ito 5/18  modified to use EXY, EXX directly to compute depm, thus
%             avoiding the interpolation of EXX to shear nodes
%=========================================================================

%dep_s=zeros(Ny,Nx);  %! this wasn't here BUT should have been! G.Ito 5/4/18

% if (PARAMS.welast==0);  %this option uses the full strain-rate 
%     dep_s(s_nodes_yield) = dt_m.*epsII_s(s_nodes_yield);
%
% else   %this option uses only the non-elastic strain
%     dep_s(s_nodes_yield) = dt_m.*max(epsII_s(s_nodes_yield)-...
%     PARAMS.welast.*(yield_s(s_nodes_yield)-sqrt(sxxOLD_s(s_nodes_yield).^2+sxyOLD(s_nodes_yield).^2))./(2.*Gs(s_nodes_yield).*dt_m),min(epsII_s(:))*1e-6);
% end
% 
% 
% [depm]=SiStER_interp_shear_nodes_to_markers(dep_s,x,y,xm,ym,icn,jcn);

%This avoids the interpolation from normal to shear nodes when computing 
% epsII_s (used in the old method above) in flow_solve. 

dEXY=0.*EXY; dEXX=EXY;  
dEXY(s_nodes_yield)=EXY(s_nodes_yield).*dt_m;
if (isfield(MAT,'psi'));
    dEXX(n_nodes_yield)=(EXX(n_nodes_yield)-0.5.*Rn(n_nodes_yield)).*dt_m;
else
    dEXX(n_nodes_yield)=EXX(n_nodes_yield).*dt_m;
end

% These are strain changes,but using (unneeded) arrays for strainrate to 
% save RAM (Exxm, EXYm are used in SiStER_update_epsIIm)
EXYm=SiStER_interp_shear_nodes_to_markers(dEXY,x,y,xm,ym,icn,jcn);  
EXXm=SiStER_interp_normal_nodes_to_markers(dEXX,xc,yc,xm,ym,icn,jcn);

dum=sqrt(EXXm.^2+EXYm.^2);

ep=(ep+dum)./(dt_m/PARAMS.tau_heal+1);
