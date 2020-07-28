%=========================================================================
% Get material properties on nodes from advected properties of markers
% G.Ito 8/16
%=========================================================================
% PHASE PROPORTIONS AT NORMAL AND SHEAR NODES. G.Ito 8/16
[n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,im);
phase_n=n2interp(1).data;
phase_n=round(phase_n*1e10)/1e10;  %prevents a case in which phase_n>NPhase

[n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,im);
phase_s=n2interp(1).data;
phase_s=round(phase_s*1e10)/1e10; %prevents a case in which phase_n>NPhase

% GET MARKER DENSITIES 
[rhom]=SiStER_get_density(im,Tm,MAT);
% pass density to nodes
[n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,rhom);
rho  = n2interp(1).data;

% GET MARKER ELASTIC PROPERTIES  G.Ito 8/16
%     [Gm]=SiStER_get_elastic_moduli(im,MAT);
%  pass shear modulus to nodes
%     [n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,1./Gm);
%     Gn=1./(n2interp(1).data);
%     [n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,1./Gm);
%     Gs = 1./(n2interp(1).data);  

[Gs]=(SiStER_get_mean_material_prop(1./cell2mat({MAT([1:Nphase]).G}),phase_s)).^-1;
[Gn(2:end,2:end)]=(SiStER_get_mean_material_prop(1./cell2mat({MAT([1:Nphase]).G}),phase_n(2:end,2:end))).^-1;

%PROPERTIES FOR PLASTICITY  G.Ito 8/16
[cohes,psim]=SiStER_get_cohesion_psim(im,ep,MAT);
[n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,cohes);
Cohes_n=n2interp(1).data;
[n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,cohes);
Cohes_s = n2interp(1).data;  

%DILATION FOR PLASTICITY
if ((PARAMS.YNPlas==1) && isfield(MAT,'psi')~=0);
    [n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,psim);
    psi_n  = n2interp(1).data;
    [n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,psim);
    psi_s  = n2interp(1).data;  %this is needed only for Rs & epsII_s in flow solve
    
    %plastic strain is needed for calculation of R
    [n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,ep);
    ep_n=n2interp(1).data;
    [n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,ep);
    ep_s = n2interp(1).data;  
end

%OPTION FOR STRAIN-DEPENDENT FRICTION
if (isfield(MAT, 'mumin')==0)
    [Mu_s]=SiStER_get_mean_material_prop(cell2mat({MAT([1:Nphase]).mu}),phase_s);
    [Mu_n(2:end,2:end)]=SiStER_get_mean_material_prop(cell2mat({MAT([1:Nphase]).mu}),phase_n(2:end,2:end));
else
    if (t==start_step)
        disp('>>>>>>>>Strain weakening for Friction<<<<<<<<')
    end
    [dum]=SiStER_get_mu(im,ep,MAT);
    [n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,dum);
    Mu_n=n2interp(1).data;
    [n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,dum);
    Mu_s = n2interp(1).data;  
end
    

if (isfield(PARAMS,'decollement_ytop'))
    SiStER_set_decollement_props;
end

%ADVECTED strainrate invariant G.Ito 8/16
[n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,epsIIm);
epsII_n=n2interp(1).data;
[n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,epsIIm);
epsII_s = n2interp(1).data;  

% OLD STRESSES AND PRESSURES G.Ito 8/16
[n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,sxxm);
sxxOLD(:,:) = n2interp(1).data;
[sxxOLD_s]=SiStER_interp_normal_to_shear_nodes(sxxOLD,dx,dy);

[n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,sxym);
sxyOLD(:,:) = n2interp(1).data;  
[sxyOLD_n]=SiStER_interp_shear_to_normal_nodes(sxyOLD);

%Now advecting pressure too!  G.Ito 5/18
[n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,pm);
pold= n2interp(1).data;
[n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,pm);
ps_old= n2interp(1).data;

% EXYOLD=EXY;
% EXXOLD=EXX;
% EXX_sOLD=SiStER_interp_normal_to_shear_nodes(EXX,dx,dy);;
% EXY_nOLD=SiStER_interp_shear_to_normal_nodes(EXY);

%TEMPERATURE ARRAYS NEEDE FOR DUCTILE RHEOLOGY  G.Ito 8/16
if (PARAMS.Tsolve==1);
    Ts=SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,Tm);
    Tn=SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,Tm);
end