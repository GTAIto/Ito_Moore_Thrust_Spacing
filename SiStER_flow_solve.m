%==========================================================================
% SiStER_flow_solve
% Performs inner solve of linear LS=R system as well as outer, iterative
% solution for non-linear dependence of L (viscosity) on S (vx,vy,P)
% Used to be named "run_Picard_iterations" but name changed by G.Ito 6/21/16
%==========================================================================
ResL2=1; 

for pit=1:maxPicard

    %% ---------------------------------------------------------------------------------
    % Compute visco-elasto-plastic viscosities
    %---------------------------------------------------------------------------------
    SiStER_VEP_rheology   %<<<<<<Interpolation p to ps done in here

    %---------------------------------------------------------------------------------
    % Assemble L and R matrices
    %---------------------------------------------------------------------------------
    [L, R, Kc, Kb]=SiStER_assemble_L_R(dx,dy,Zs.*etas,Zn.*etan,rho,BC,PARAMS,srhs_xx,srhs_xy,Rn,psurf); %G.Ito
    
    %---------------------------------------------------------------------------------
    % Residual:  L and R are from current solution S
    %---------------------------------------------------------------------------------
    if (exist('S','var'));
        Res=L*S-R;
        ResL2=norm(Res,2)/norm(R,2);
    end;
    %---------------------------------------------------------------------------------
    % Solve for new solution S using Picard or approximate Newton or a
    % combination of the two 
    %---------------------------------------------------------------------------------  
    %if(pit>2 && ResL2<0.1);
    if(pit > 1e6);
       beta=1;
       S=S-beta.*(L\Res);     %approximate Newton update, with L and approximation to J
       it_type='Newton: ';
    else;
       S=L\R;       %Picard update
       it_type='Picard: ';
    end;
    
    if (max(abs(imag(S))) >0)  %debugging
       disp('S is imaginary');
       crashed
    end;
    
    [p, vx, vy]=SiStER_reshape_solver_output(S,Kc,Nx,Ny);
    
    %% ---------------------------------------------------------------------------------
    % new strainrates
    % This definition of epsIIm recovers the Mohr-Coloumb yield stress
    % Rn is one iteration behind here. 
    %---------------------------------------------------------------------------------

    if (isfield(BC,'bot_xbackstop') && BC.bot_xbackstop(2)>0);
        [EXX,EXY]=SiStER_get_strain_rate_backstop(vx,vy,dx,dy,x,BC);
    else
        [EXX,EXY]=SiStER_get_strain_rate(vx,vy,dx,dy,BC);
    end;
    
    EXY_n=SiStER_interp_shear_to_normal_nodes(EXY);       %<<<< Interpolation done here
    EXX_s=SiStER_interp_normal_to_shear_nodes(EXX,dx,dy); %<<<< Interpolation done here
    if (isfield(MAT,'psi'));
        epsII_n=sqrt(0.5*EXX.^2  + 0.5*(0.5*Rn-EXX).^2  + EXY_n.^2 );  
        epsII_s=sqrt(0.5*EXX_s.^2 +0.5*(0.5*Rs-EXX_s).^2 + EXY.^2 );
    else
        epsII_n=sqrt(EXX.^2  + EXY_n.^2);  
        epsII_s=sqrt(EXX_s.^2 + EXY.^2);
    end

    %% ---------------------------------------------------------------------
    % TO EVALUATE CONVERGENCE
    %---------------------------------------------------------------------
           
    % get velocity field from current solution on shear nodes
    SiStER_interp_velocities_to_shear_nodes;  
    %dv=abs(sqrt(vxc(:).^2+vyc(:).^2)-v(:))./median(v(:));       
    dvnorm=norm(abs(sqrt(vxc(:).^2+vyc(:).^2)-v(:)),2)./norm(v);
    v=sqrt(vxc.^2+vyc.^2); 
    % fraction of domain size where change in velocity field is greater than
    
    %----------------------------------------------------------------------
    % ASSESS CONVERGENCE HERE
    %----------------------------------------------------------------------
    if (t == 1 && pit==1)
        outit=[0 0 0 1 0];
        iteration_count=1;
        conv_rate=1e6;
        nfit=2;  %sometimes ResL2 jumps between individual iterations
    else
        if (iteration_count >= nfit);
            junk=polyfit([1:nfit]',outit(end-nfit+1:end,1),1);
            conv_rate=abs(junk(1))./mean(outit(end-nfit+1:end,1));
        else
            conv_rate=1e6;
        end;
    
        outit=[outit; [ResL2 dvnorm conv_rate t pit];];
        iteration_count=iteration_count+1;
    end;
    
    picstring=([it_type sprintf('%0.f %0.4e %0.4e %0.4e %0.f %0.2f',[pit ResL2 conv_rate dvnorm t time/(3600*24*365.25*1e3)])]);
    disp(picstring);

    if(ResL2<PARAMS.conv_crit_ResL2 && pit >= PARAMS.Npicard_min)
        disp(['    Converged to conv_crit_ResL2: ' picstring]);
        break;
    elseif (conv_rate < PARAMS.conv_crit_rate && pit >= PARAMS.Npicard_min)
        disp(['Converged ONLY by conv_crit_rate: ' picstring]);
        break;
	elseif (pit==PARAMS.Npicard);
       	%disp([num2str(pit) ' Picard iterations failed to converge.']);
        disp(['               **NOT** Converged: ' picstring]);
    end
    
	v=sqrt(vxc.^2+vyc.^2);    
    
    
end

      
