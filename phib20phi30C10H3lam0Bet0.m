% SiStER_Input_File
clear all;

restart = 0;            %
restart_timestep = 800;   %only used when restart=1 
if (restart ==1)
     eval(['load restart_' num2str(restart_timestep) '.mat']); 
     restart = 1;
     start_step=t+1;
     [Ny Nx] = size(p);   
     disp(['****Restarting on time step: ' num2str(t) ' ****']) 
end;


%% DURATION OF SIMULATION AND FREQUENCY OF OUTPUT 
Nt     = 10; % max number of time iterations
dt_out = 10; % output files every "dt_out" iterations
dt_out_restart = 200;     % number of steps between outputing everything.  

isee_timesteps = 5;  % = 0 to not view output realtime, if >0
                     %then view output every isee_timesteps
See_timestep_script='Visualize_Output';
                     
                     
%% DOMAIN SIZE AND GRIDDING 
%GRID.x defines the points in x within which grid spacing can differ
%for x<=GRID.x(1), the APPROXIMATE grid spacing is GRID.dx(1); 
%for GRID.x(1) < x <= GRID.x(2), the approx. spacing is GRID.dx(2);
%for GRID.x(end-1) < x <= xsize, the approx. spacing is GRID.dx(end);
%
%Thus GRID.dx MUST be longer than GRID.x by 1 entry.
%GRID.x=[0] and GRID.dx=[0.5 1]*1e3 will sent a UNIFORM spacing of 1000m
%The grid is defined in SiStER_initialize_grid_GI.  This is also 
%where Nx will be defined.
%
%To ensure an EXACT grid spacing, be sure GRID.dx fits as an integer number
%of times in each grid interval GRID.x

%The same logic is used for y.
%--------------------------------------------------------------------

xsize = 65.e3;
ysize = 11.0e3;
%GRID.x=[0]; GRID.dx=[1 1]*1e3;  %Uniform grid spacing at 1000m
%GRID.y=[0]; GRID.dy=[1 1]*1e3;  %Uniform grid spacing at 1000m
GRID.x=[40.0]*1e3;                  %gridding domains
GRID.dx=[0.1 0.1]*1e3;         %make sure these fit an integer number of
                                 %times in each LAYER
GRID.y=[10.8]*1e3;                  
GRID.dy=[0.1 0.1]*1e3;

% LAGRANGIAN MARKERS 
Mquad      = 8; % number of markers in the smallest quadrant
Mquad_crit = 2; % minimum number of markers allowed in smallest quadrant (for reseeding)

%% GEOMETRY OF MATERIALS

Nphase=2; % number of phases

% phase 1
GEOM(1).type = 1; % 1 = layer (then specify top and bot) or 2 = circle 
GEOM(1).top  = 0;
GEOM(1).bot  = 8.0e3;

% phase 2
GEOM(2).type = 1; % 1 = layer (then specify top and bot) or 2 = circle 
GEOM(2).top  = GEOM(1).bot;
GEOM(2).bot  = 11.0e3;


% PASSIVE STRATIGRAPHIC LAYER LOCATIONS FOR VISUALIZING DEFORMATION
strata_y=[GEOM(2).top:0.2e3:GEOM(2).bot];

%% MATERIAL PROPERTIES 

% phase 1
MAT(1).phase=1;
% density parameters
MAT(1).rho0=1000.0;
MAT(1).alpha=0;
% elasticity 
MAT(1).G=5e10;
% diffusion creep parameters
MAT(1).pre_diff=1/1e17;
MAT(1).Ediff=0;
MAT(1).ndiff=1;
MAT(1).mdiff=0;
% dislocation creep parameters
MAT(1).pre_disc=1/1e17;
MAT(1).Edisc=0;
MAT(1).ndisc=1;
MAT(1).mdisc=0;
% grain size
MAT(1).dgrain=0;
% plasticity
MAT(1).mu=tand(30);
MAT(1).Cmax=10.0e6;
MAT(1).Cmin=1e6;
MAT(1).ep_startweakening=0.002;
MAT(1).ecrit=0.1;
% thermal
MAT(1).k=0.001;
MAT(1).cp=0.1;


% phase 2
MAT(2).phase=2;
% density parameters
MAT(2).rho0=2300;
MAT(2).alpha=0;
% elasticity 
MAT(2).G=5e10;
% diffusion creep parameters
MAT(2).pre_diff=1/1e24;
MAT(2).Ediff=0;
MAT(2).ndiff=1;
MAT(2).mdiff=0;
% dislocation creep parameters
MAT(2).pre_disc=1/1e35;
MAT(2).Edisc=0;
MAT(2).ndisc=1;
MAT(2).mdisc=0;
% grain size
MAT(2).dgrain=0;
% plasticity
MAT(2).mu=tand(30);
MAT(2).Cmax=10.0e6;
MAT(2).Cmin=1e6;
MAT(2).ep_startweakening=0.002;
MAT(2).ecrit=0.1;
MAT(2).psi=atand(MAT(2).mu);
MAT(1).psi=atand(MAT(2).mu);
% thermal
MAT(2).k=0.001;
MAT(2).cp=0.1;

PARAMS.decollement_ytop = GRID.y(1);
PARAMS.decollement_mu=tand(20);
PARAMS.decollement_coh=0.0;

%% ADDITIONAL PARAMETERS 
PARAMS.YNElast = 1; % elasticity on (1) or off (0)
PARAMS.YNPlas  = 1; % plasticity on (1) or off (0)

PARAMS.tau_heal     = 0.4e12; % healing time for plasticity (s)
% PARAMS.plastic_seed = 'random'; % use plastic strain seed to promote faulting
% PARAMS.plastic_seed = 'box   ';
% Seed initial weak zone(s)
% if (PARAMS.plastic_seed=='box   ');
%     PARAMS.seed_dim=750;
%     PARAMS.seed_x=5e3;
%     PARAMS.seed_y=8e3;
%     PARAMS.seed_amp=0;
% end;

beta=0.0;   %Basal dip
PARAMS.gx = 9.8*sind(beta); % gravity along x
PARAMS.gy = 9.8*cosd(beta); % gravity along y
PARAMS.ridge=0;
PARAMS.fracCFL    = 0.5;     % distance by which a marker is allowed to move 
                             % over a time step, as a fraction of the smallest cell size
PARAMS.grainscale = 1;       % to normalize the grain size in the creep laws
PARAMS.R          = 8.314;   % gas constant
PARAMS.etamax     = 1e24;    % maximum viscosity
PARAMS.etamin     = 1e17;    % minimum viscosity
PARAMS.plast_option = 2; % keep full strainrate (dont remove elastic part)
PARAMS.Npicard    = 5;      % number of Picard iterations per time step
PARAMS.Npicard_min    = 5;      % number of Picard iterations per time step
PARAMS.conv_crit_ResL2 = 0.001;   %non-linear residual norm defining convergence G.Ito     
PARAMS.conv_crit_rate = 0.0002;  %min rate of change in ResL2 defining convergence G.Ito
PARAMS.rem_water_col_press = 1;  %remove surface water column pressure
PARAMS.lambda_pp=0.0;  %pore-fluid pressure ratio

PARAMS.p0cell     = 0;       % Pressure in pressure anchor
PARAMS.Tsolve     = 0;       % yes (1) or no (0) solve for temperature

% TOPOGRAPHY
isurface  = 1;                %surface topography (x_topo, y_topo) now marked
PARAMS.dxsurf     = min(GRID.dx)./4; 
PARAMS.kappa      = 6.25e-9;    %surface topography diffusion
PARAMS.reseedsurface = 2;    %reseed surface markers between this many timesteps

% initial temperature profile, polynomial with depth 
% T = a0 + a1*y+a2*y^2+a3*y^3+amp*sin(2*pi*X/lam)
% (make sure it matched the BCs)
PARAMS.a0  = 0;
PARAMS.a1  = 0;
PARAMS.a2  = 0;
PARAMS.a3  = 1;
PARAMS.amp = 0;
PARAMS.lam = 1;

%--------------------------------------------------------------------------------------
%% BOUNDARY CONDITIONS 

% Flow boundary conditions
% entries in BC correspond to [a b c d]  %G.Ito added "d"
% a/ rollers? 
%    0 = No slip
%    1 = Free slip (if using multiple side BCs, set in BC.side, not
%        BC.side_w/l/m
% b/ type of velocity normal to boundary 
%    0 = constant, cosnervative 
%    1 = triangular profile (top/bot only)
%    2 = multiple (sides only)
%    3 = open     (only for top, bottom, and multiple)
% c/ value of normal velocity (or max. value along the profile, if b/ is 1)
%
% d/ value of tangential velocity for No Slip (only works for bottom boundary)

Uin = 1e-2/(365.25*24*3600);

BC.top     = [1  3  0];      
BC.bot     = [0  0  0 -Uin];  %4th number is tangential vel. for no slip
%BC.bot_xbackstop = [5.0e3 5.5e3];  % backstop x < (1) U=0, x>(2) U=BC.bot(4);

BC.left    = [0 0    0];
BC.right   = [1 0  -Uin]; 

% Thermal 
% ------------------
% entries in BCtherm correspond to
% 1/ type? (1=Dirichlet, 0=Neumann)
% 2/ value
BCtherm.top   = [1 0];
BCtherm.bot   = [1 1];
BCtherm.left  = [2 0];
BCtherm.right = [2 0];


%--------------------------------------------------
% Run it
%--------------------------------------------------
SiStER_MAIN;

