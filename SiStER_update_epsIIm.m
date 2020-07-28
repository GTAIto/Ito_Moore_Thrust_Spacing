%----------------------------------------------------------------------------
%Put strain-rate invariant on markers so they will be advected.
% G.Ito 5/7/18 Updated to use EXX and EXY directly rather than using epsII_s (which
% involved interpolating EXX to shear nodes in flow solve). 
% This definition of epsIIm recovers the Mohr-Coloumb yield stress
%---------------------------------------------------------------------------

%   Old code
%    epsIIm=SiStER_interp_shear_nodes_to_markers(epsII_s,x,y,xm,ym,icn,jcn);

%New code
 EXYm=SiStER_interp_shear_nodes_to_markers(EXY,x,y,xm,ym,icn,jcn);
 if (isfield(MAT,'psi')~=0)
    EXXm=SiStER_interp_normal_nodes_to_markers((EXX-0.5.*Rn),xc,yc,xm,ym,icn,jcn);
 else
    EXXm=SiStER_interp_normal_nodes_to_markers(EXX,xc,yc,xm,ym,icn,jcn);
 end
 epsIIm=sqrt(EXXm.^2+EXYm.^2); 