% SiStER THERMAL SOLVE

% get previous temperature on nodes
[n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,Tm);
Told(:,:)  = n2interp(1).data;


% enforce Dirichlet boundary conditions to avoid mismatch between markers
% and nodes
if BCtherm.top(1)==1
    Told(1,:)=BCtherm.top(2);
end
if BCtherm.bot(1)==1
    Told(Ny,:)=BCtherm.bot(2);
end
if BCtherm.left(1)==1
    Told(:,1)=BCtherm.left(2);
end
if BCtherm.right(1)==1
    Told(:,Nx)=BCtherm.right(2);
end

% get thermal conductivity and heat capacity on markers
[km, cpm]=SiStER_get_thermal_properties(im,MAT);
% we already have rho on shear nodes
% pass heat capacity and thermal conductivity to shear nodes
[n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,qd,x,y,cpm,km);
cp(:,:)  = n2interp(1).data;
k(:,:)   = n2interp(2).data;
% pass thermal conductivity to vy (kx) and vx (ky) nodes   
kx=zeros(size(k));
ky=zeros(size(k));
kx(:,1:end-1)=0.5*(k(:,1:end-1)+k(:,2:end));
ky(1:end-1,:)=0.5*(k(1:end-1,:)+k(2:end,:));

[T]=SiStER_thermal_solver_sparse(x,y,Told,rho,cp,kx,ky,dt_m,BCtherm);  
% temperature change
dT=T-Told;
% enforce Dirichlet boundary conditions to avoid mismatch between markers
% and nodes
if BCtherm.top(1)==1
    dT(1,:)=0;
end
if BCtherm.bot(1)==1
    dT(Ny,:)=0;
end
if BCtherm.left(1)==1
    dT(:,1)=0;
end
if BCtherm.right(1)==1
    dT(:,Nx)=0;
end

% pass temperature change to markers
[dTm]=SiStER_interp_shear_nodes_to_markers(dT,x,y,xm,ym,icn,jcn);
% update marker temperatures
Tm=Tm+dTm;