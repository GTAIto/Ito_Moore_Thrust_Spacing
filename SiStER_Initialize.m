%===========================================================================
% Create arrays and initalize their values
%===========================================================================


% construct staggered grids
[X,Y,x,y,xc,yc,dx,dy,Nx,Ny] = SiStER_initialize_grid_GI(xsize,ysize,GRID);  %G.Ito
    
% initialize marker arrays and positions
[xm, ym] = SiStER_initialize_marker_positions(xsize,ysize,dx,dy,Mquad);

% locate markers with respect to grid
[qd,icn,jcn] = SiStER_locate_markers_in_grid(xm,ym,x,y,dx,dy);

% assign marker phases
if (isfield(BC,'bot_xbackstop') && BC.bot_xbackstop(2)>0);
    %[im] = SiStER_initialize_marker_phases_backstop(Nphase,GEOM,xm,ym,PARAMS,xsize);
    disp('BACK_STOP not working in this version')
    halt
else
    [im] = SiStER_initialize_marker_phases(Nphase,GEOM,xm,ym,PARAMS,xsize);
end
%% INITIALIZE ARRAYS
% initialize marker plastic strain (to zero) and strain rate (to one)
ep=zeros(size(xm));
epsIIm=ones(size(xm))*1e-20;

% initialize marker stresses
sxxm=zeros(size(xm));
sxym=sxxm;
etam=sxxm;
dsxxm=sxxm;
dsxym=sxym;
pm=sxxm;    %Added G.Ito 5/18
dpm=sxxm;   %Added G.Ito 5/18

% initialize marker index (a unique number to identify and track each marker)
idm=1:length(xm);

% initialize stratigraphic layer numbers (passive markers of deformation) G.Ito
if (exist('strata_y','var'))                                        %G.Ito
    [stratam]=SiStER_initialize_marker_stratigraphy(strata_y,xm,ym); %G.Ito
else
    stratam=zeros(size(sxxm));                                      %G.Ito
end                                                                %G.Ito

% initialize temperature stucture on nodes
if (PARAMS.Tsolve)
    T=PARAMS.a0+PARAMS.a1*Y+PARAMS.a2*Y.^2+PARAMS.a3*Y.^3;
    T=T+PARAMS.amp*sin(2*pi*X/PARAMS.lam);
    % pass initial nodal T to markers
    [Tm]=SiStER_interp_shear_nodes_to_markers(T,x,y,xm,ym,icn,jcn);
else
    Tm=zeros(size(xm));
    Ts=zeros(Ny,Nx);
    Tn=zeros(Ny,Nx);
end

EXX=zeros(size(X));
EXY=zeros(size(X));
epsII_s=1e-36.*ones(size(EXY));  %make it small but non-zero for eta_plas
epsII_n=1e-36.*ones(size(EXX));
dep_s=zeros(Ny,Nx);
etan_new=zeros(Ny,Nx);
yield_n=zeros(Ny,Nx);
Mu_n=zeros(Ny,Nx);
Gn=zeros(Ny,Nx);
Gs=zeros(Ny,Nx);

p=1e9*ones(size(EXX));  %initialize to be high so plasticity doesnt activate at t=1, pit=1;

% initialize velocities (for convergence test in Picard iterations) G.Ito 
vx=EXX; vy=EXX;

Rn=zeros(Ny,Nx); %Dilation w/plasticity: this must be initialized for assemble_L_R
Rs=0.0; %for output

%%
%Set minimum number of Picard iterations
if isfield(PARAMS,'Npicard_min')==0;
         PARAMS.Npicard_min=1;
end

% FAULTING SEED
% randomly seeds ep to promote faulting if switch is on (internal)
if (exist('PARAMS.plastic_seed','var') && sum(PARAMS.plastic_seed ~= 'none  ')<6);
    SiStER_plastic_seed;
end

%% -------------------------------------------------------------------------
% SURFACE MARKERS  G.Ito
%-------------------------------------------------------------------------
if exist('isurface','var')
    if isurface==1
        Ntopo=round(xsize/PARAMS.dxsurf)+1;
        topo_x=linspace(0,xsize,Ntopo);
        topo_y=GEOM(1).bot*ones(size(topo_x));
        if (isfield(BC,'bot_xbackstop') && BC.bot_xbackstop(2)>0)
            kk=4;
            dum1=GEOM(kk).ytopleft+(GEOM(kk).yright-GEOM(kk).ytopleft)/(GEOM(kk).xright-GEOM(kk).xtopleft)*(topo_x-GEOM(kk).xtopleft);
            dum1=min([dum1' GEOM(kk).yright.*ones(size(topo_y')) topo_y']');
            topo_y=dum1;
        end
    end
end

%% -------------------------------------------------------------------------
% initialize dt_m small to keep things elastic & no plasticity at t=1, G.Ito
%-------------------------------------------------------------------------
if (exist('dt_m','var')==0)
    dt_m=1e2;
end
sxxOLD=EXX; sxyOLD=EXY;  %inialize these as zeros G.Ito


%% --------------------------------------------------------------------------
% Set DEFAULT params G.Ito
%--------------------------------------------------------------------------
if (isfield(PARAMS,'plast_option')==0)
    PARAMS.plast_option=2;           %1=harmonic sum of eta_plas and eta_new, 2=limit eta_new to eta_plas in VEP_rheology
end

if (isfield(PARAMS,'welast')==0)
    PARAMS.welast=0;           %1=harmonic sum of eta_plas and eta_new, 2=limit eta_new to eta_plas in VEP_rheology
end

if (isfield(BC, 'bot_taperx')==0)
    BC.bot_taperx=[x(5) x(end-5)];   %default is 5 elements from the end
end

if (isfield(BC, 'top_taperx')==0)
    BC.top_taperx=[x(5) x(end-5)];   %default is 5 elements from the end
end

%% -------------------------------------------------------------------------
% Open top for tilted model (gx~=0) 
% psurf is needed for the top boundary if MAT(1).rho0>0;
%---------------------------------------------------------------------------
psurf=zeros(1,Nx);
if (BC.top(2)==3 && PARAMS.gx < 0)
    tilt=atand(PARAMS.gx/PARAMS.gy);   
    disp(['**** Model titled counter clockwise: tilt=' num2str(tilt) ' degrees']);
    rhog=MAT(1).rho0.*sqrt(PARAMS.gx.^2+PARAMS.gy.^2);
    psurf(2:Nx)=((X(1,1:Nx-1)+X(1,2:Nx))/2-X(1,Nx)).*sind(tilt).*rhog;
elseif (BC.top(2)==3 && PARAMS.gx > 0)
    tilt=atand(PARAMS.gx/PARAMS.gy);   
    disp(['**** Model titled clockwise: tilt=' num2str(tilt) ' degrees']);
    rhog=MAT(1).rho0.*sqrt(PARAMS.gx.^2+PARAMS.gy.^2);
    psurf(2:Nx)=((X(1,1:Nx-1)+X(1,2:Nx))/2).*sind(tilt).*rhog;
end
%% -------------------------------------------------------------------------
% Option to remove water column
% Turn on only if pore-pressure via lambda is off
if (isfield(PARAMS,'rem_water_col_press')==0)
    PARAMS.rem_water_col_press=0;
end
%% -------------------------------------------------------------------------
% Hubbert-Rubey pore pressure option
if (isfield(PARAMS,'lambda_pp')== 0)
    PARAMS.lambda_pp=0;
end
if (PARAMS.lambda_pp > 0 || PARAMS.rem_water_col_press > 0)
    dYn=repmat([0 0.5*dy(1) (dy(2:end)+dy(1:end-1))./2]',1,Nx);  %cell centers for Plith calculation in VEP_rheology
    dYs=repmat(dy',1,Nx); 
    disp(['Water removed and/or lambda on: lambda_pp, rem_water_col=' num2str([PARAMS.lambda_pp PARAMS.rem_water_col_press])]);
end

if (PARAMS.lambda_pp>0 && PARAMS.rem_water_col_press>0)
    fdum=zeros(Ny,Nx);  %initialize array for fraction of water in element
end




%% -------------------------------------------------------------------------
% Material property OPTIONS
%-------------------------------------------------------------------------
nonNewtonian=0;
for k=1:Nphase
    if (MAT(k).ndiff ~=1 || MAT(k).ndisc ~=1)
        disp('Sorry cannot do power-law rheology yet')
        halt
    end;
end
 if (isfield(MAT,'ep_startweakening')==0)
    for k=1:Nphase   
        MAT(k).ep_startweakening=0;
    end
 else
    for k=1:Nphase  
        dum=MAT(k).ep_startweakening;
        if (length(dum)==0)
            MAT(k).ep_startweakening=0;
        end
    end
 end
%% -------------------------------------------------------------------------
% Display parameters to std out
%-------------------------------------------------------------------------
for k=1:Nphase
	MAT(k)
end
    
PARAMS

BC



