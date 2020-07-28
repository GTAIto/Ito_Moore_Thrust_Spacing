% SiStER_MAIN.m
%
% Simple Stokes with Exotic Rheologies
%
% Main routine doing initialization, time loop and outputs
%
%
% J.-A. Olive, B.Z. Klein, E. Mittelstaedt, November 2014
% S.M. Howell June 2015
% G.Ito 2015-2016
% See github.com/jaolive/SiStER
% clear all
% close all

% INITIALIZATION

% % load parameter values, model geometry, boundary conditions
% input_file  = 'SiStER_Input_File_extension';  %G.Ito
% output_name = 'test_3pic';  %G.Ito

% Load it
% eval(input_file);
tic;
c = clock;
% if (exist('output_dir','var')==0);
%     error('output_dir not specified');
% end;
output_dir = '.';
if restart==0                                                %G.Ito/S.Howell
    SiStER_Initialize
    time=0;
    start_step = 1;                                         %G.Ito/S.Howell
end

% For TANGENTIAL VELOCITY BCs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SiStER_BC_velocity_profiles;

% BEGIN TIME LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t=start_step:Nt

    time=time+dt_m;  %Time of THIS step is with current temperatures and dt_m used in VEP_rheology
    tcpu0=cputime; dcpu_markers=0; dcpu_stokes=0; %G.Ito 

    SiStER_material_props_on_nodes;


    %% --------------------------------------------------------------------------
    % PRODUCE CURRENT SOLUTION FOR FLOW, STRESS, and STRAINRATE
    % --------------------------------------------------------------------------
    tcpu1=cputime;  dcpu_markers=tcpu1-tcpu0; tcpu_stokes=0; %G.Ito
    if t==1
        v=vx;          %just to to check convergence in Picard iteration
        maxPicard=PARAMS.Npicard_min;
    else
        maxPicard=PARAMS.Npicard;
    end
    if (PARAMS.YNPlas==0 && nonNewtonian==0);
        maxPicard=1;
    end

    SiStER_flow_solve;

    % UPDATE MARKER STRESSES FOR CURRENT SOLUTION (ROTATION OCCURS BELOW)
    if (PARAMS.YNElast==1) 
        SiStER_update_marker_stresses;
    end
    tcpu2=cputime;  %G.Ito  

    %----------------------------------------------------------------------------
    %Put strain-rate invariant on markers so they will be advected.
    % G.Ito 5/7/18 Updated to use EXX and EXY directly rather than using epsII_s (which
    % involved interpolating EXX to shear nodes in flow solve). 
    % This definition of epsIIm recovers the Mohr-Coloumb yield stress
    %----------------------------------------------------------------------------
    if (PARAMS.YNPlas==1 || nonNewtonian==1)
        SiStER_update_epsIIm;
    end;

    %% -------------------------------------------------------------------------- 
    % OUTPUT solution of current time step
    %-------------------------------------------------------------------------- 
    SiStER_output;
    % >>>REAL TIME VISUALIZATION G.Ito
    if (exist('isee_timesteps','var') && (mod(t,isee_timesteps)==0));
        eval(See_timestep_script);
        disp('Hit return to continue');
        pause;
    end

    %% ------------------------------------------------------------------------
    % Now evolve properties to the next timestep

    % SET TIME STEP for updates to stresses, temperature, and markers
    % needed for next timestep etc. G.Ito
    if (exist('dt_m_static','var')==1);
        dt_m=dt_m_static;
    else
        [dt_m]=SiStER_set_timestep(dx,dy,vx,vy,PARAMS);
    end
    % EVOLVE ACCUMULATED PLASTIC STRAIN
    if (PARAMS.YNPlas==1) %G.Ito
        SiStER_update_ep;  
    end

    % ROTATE ELASTIC STRESSES for NEXT STEP
    if (PARAMS.YNElast==1) 
        SiStER_rotate_stresses;
    end


    dcpu_markers=dcpu_markers+cputime-tcpu2;

    %--------------------------------------------------------------------------

    % THERMAL SOLVE  
    if PARAMS.Tsolve==1
        SiStER_thermal_update;
    end


    tcpu3=cputime;  %G.Ito
    % MARKER UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SiStER_update_markers;
    % advect markers
    % remove markers if necessary
    % add markers if necessary
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %>>epsIIm now passed through SiStER_update_markers
    %get current nodal strain rate on updated markers
    %  [exxm]=SiStER_interp_normal_nodes_to_markers(EXX,xc,yc,xm,ym,icn,jcn);
    %  [exym]=SiStER_interp_shear_nodes_to_markers(EXY,x,y,xm,ym,icn,jcn);
    %  epsIIm=sqrt(exxm.^2+exym.^2);

    dcpu_markers=cputime-tcpu3+dcpu_markers;  %G.Ito

    %Surface topography markers and erosion
    tcpu4=cputime; %G.Ito
    SiStER_update_surface; %G.Ito

    dcpu_surface=cputime-tcpu4;

    % end time loop, update time


    disp('===================================')
    disp(['END TIMESTEP=' num2str(t) '; time=' num2str(time/(365.25*24*60*60*1e6)) ....
        ' (Myrs); # Picards=' num2str(pit),'; Errors= ', num2str([ResL2 dvnorm])]);
    disp(['*CPU total, markers, Stokes, non-Linear iterations, Surface=' ....
          num2str((cputime-tcpu0)/60) ', ' num2str(dcpu_markers/60) ', ' ....
          num2str(dcpu_stokes/60) ', ' num2str((tcpu2-tcpu1)/60) ', ' num2str(dcpu_surface/60)]);  %G.Ito
    disp('===================================')

    %if (mod(t,5)==0); SiStER_see; disp('Hit any key to continue'); pause; end;
    if (exist('isurface','var')==1 && isurface==1)
    if (min(topo_y)<y(2))
        disp('STOPPING: Wedge Surface at top of box')
        break
    end
    end

end

toc;

    
