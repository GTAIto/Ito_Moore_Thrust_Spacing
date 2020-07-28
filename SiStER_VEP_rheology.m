%==========================================================================
% SiStER_UPDATE_RHEOLOGY
% calculates visco-elasto-plastic rheology terms on shear and normal nodes
% G.Ito 8/20
%==========================================================================

%--------------------------------------------------------------------------
% plastic viscosity first (Mohr-Coulomb)
%--------------------------------------------------------------------------
if (pit==1)
    ps=ps_old;  %ps_old is from pm, just advected from previous t-step G.Ito 5/18
    if ((PARAMS.lambda_pp>0 || PARAMS.rem_water_col_press==1) && PARAMS.YNPlas==1)  %For pore pressure or water column pressure
        Plithn=[zeros(1,Nx); [zeros(Ny-1,1) (rho(1:Ny-1,1:Nx-1)+rho(1:Ny-1,2:Nx))./2];];  %densities above normal & shear nodes
        Pliths=[zeros(1,Nx); (rho(1:Ny-1,:)+rho(2:Ny,:))./2;];
        
        if (PARAMS.rem_water_col_press==1 && PARAMS.lambda_pp==0)  
            %(water fraction)*(water density)
            Plithn(Plithn<MAT(2).rho0)=MAT(1).rho0.*(MAT(2).rho0-Plithn(Plithn<MAT(2).rho0))./(MAT(2).rho0-MAT(1).rho0);  
            Plithn(Plithn>=MAT(2).rho0)=0;  %no rock weight to remove
            Pliths(Pliths<MAT(2).rho0)=MAT(1).rho0.*(MAT(2).rho0-Pliths(Pliths<MAT(2).rho0))./(MAT(2).rho0-MAT(1).rho0);
            Pliths(Pliths>=MAT(2).rho0)=0;
        elseif(PARAMS.rem_water_col_press==1 && PARAMS.lambda_pp>0)   
            %divide by PARAMS.lambda_pp to keep 100% of the water column
            fdum=(MAT(2).rho0-Plithn)./(MAT(2).rho0-MAT(1).rho0); fdum(fdum<0)=0; fdum(fdum>1)=1;  %fraction water layer
            Plithn(Plithn<MAT(2).rho0)=MAT(1).rho0.*fdum(Plithn<MAT(2).rho0)./PARAMS.lambda_pp+... %water layer
                                       MAT(2).rho0.*(1-fdum(Plithn<MAT(2).rho0));                  %rock
                                   
            fdum=(MAT(2).rho0-Pliths)./(MAT(2).rho0-MAT(1).rho0); fdum(fdum<0)=0; fdum(fdum>1)=1;
            Pliths(Pliths<MAT(2).rho0)=MAT(1).rho0.*fdum(Pliths<MAT(2).rho0)./PARAMS.lambda_pp+...
                                       MAT(2).rho0.*(1-fdum(Pliths<MAT(2).rho0));
            
        end
        Plithn(2:Ny,:)=cumsum(Plithn(2:Ny,:).*dYn(2:Ny,:)).*PARAMS.gy;        
        Pliths(2:Ny,:)=cumsum(Pliths(2:Ny,:).*dYs).*PARAMS.gy;
    end

else
    dp=p-pold;  %update ps from previous Picard iteration
    ps=ps_old+SiStER_interp_normal_to_shear_nodes(dp,dx,dy);  %<<<< Interpolation done here
end

if (PARAMS.YNPlas==1)
    pnn=p; 
    if (PARAMS.lambda_pp>0)   %Hubbert-Rubey pore pressure
         ps=ps-PARAMS.lambda_pp.*Pliths;
         pnn=pnn-PARAMS.lambda_pp.*Plithn;
    elseif (PARAMS.rem_water_col_press==1)
        pnn=pnn-Plithn;
        ps=ps-Pliths;
    end
    ps(ps<0)=0;
    pnn(pnn<0)=0;
    yield_s=(Cohes_s+Mu_s.*ps).*cos(atan(Mu_s));
    yield_n=(Cohes_n+Mu_n.*pnn).*cos(atan(Mu_n));
    if (PARAMS.welast==0)
        eta_plas_s=0.5.*yield_s./epsII_s;  
        eta_plas_n=0.5.*yield_n./epsII_n;
    else
        eta_plas_s=0.5.*yield_s./max(epsII_s-PARAMS.welast*(yield_s-sqrt(sxxOLD_s.^2+sxyOLD.^2))./(2.*Gs.*dt_m),min(epsII_s(:))*1e-6);  
        eta_plas_n=0.5.*yield_n./max(epsII_n-PARAMS.welast*(yield_n-sqrt(sxxOLD.^2+sxyOLD_n.^2))./(2.*Gn.*dt_m),min(epsII_n(:))*1e-6); 
    end
end

%-------------------------------------------------------------------------
% Ductile viscosity
%-------------------------------------------------------------------------
[etas_new]=SiStER_get_ductile_rheology(MAT,PARAMS,Ts,epsII_s,phase_s);
[etan_new(2:end,2:end)]=SiStER_get_ductile_rheology(MAT,PARAMS,Tn(2:end,2:end),epsII_n(2:end,2:end),phase_n(2:end,2:end));

if (PARAMS.YNPlas==1);
    s_nodes_yield=find(etas_new>eta_plas_s);
    n_nodes_yield=find(etan_new>eta_plas_n);
    if (PARAMS.plast_option==1);  
        etan_new(2:end,2:end)=(1./eta_plas_n(2:end,2:end) + 1./etan_new(2:end,2:end) + 1./PARAMS.etamax).^-1;
        etas_new=(1./eta_plas_s + 1./etas_new + 1./PARAMS.etamax).^-1;
        if (pit==1); disp('Using harmonic mean for eta_vp'); end; 
    else
        etan_new(2:end,2:end)=(1./etan_new(2:end,2:end) + 1./PARAMS.etamax).^-1;
        etas_new=(1./etas_new + 1./PARAMS.etamax).^-1;

        etan_new(n_nodes_yield)=eta_plas_n(n_nodes_yield);
        etas_new(s_nodes_yield)=eta_plas_s(s_nodes_yield);
        if (pit==1); disp('Using sharp plasticity cut-off for eta_vp'); end;
    end
    %if(isfield(MAT,'psi') && pit==1);
    if(isfield(MAT,'psi'));
        Rn=zeros(Ny,Nx); Rs=Rn;
        Rn(n_nodes_yield)=2.*epsII_n(n_nodes_yield).*sind(psi_n(n_nodes_yield));
        Rs(s_nodes_yield)=2.*epsII_s(s_nodes_yield).*sind(psi_s(s_nodes_yield));
    end;
end


alpha=1;  %weighting of update to minimize nodes switching in/out of failure (does not seem to be needed)
if (pit==1);
    etan=etan_new;
    etas=etas_new;
else
    etan=(etan_new.^alpha).*(etan.^(1-alpha));
    etas=(etas_new.^alpha).*(etas.^(1-alpha));
end
etan(etan<PARAMS.etamin)=PARAMS.etamin;
etas(etas<PARAMS.etamin)=PARAMS.etamin;
%-------------------------------------------------------------------------
% ELASTICITY TERMS 
%-------------------------------------------------------------------------
if PARAMS.YNElast==0
    Zs=ones(size(etas));
    Zn=ones(size(etan));
else
    Zs=Gs*dt_m./(Gs*dt_m+etas);
    Zn=Gn*dt_m./(Gn*dt_m+etan);
end
% right-hand size (1-Z)*sigmaOLD_ij
srhs_xx=(1-Zn).*sxxOLD;
srhs_xy=(1-Zs).*sxyOLD;



