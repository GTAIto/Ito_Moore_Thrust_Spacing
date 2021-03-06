erosion_by_diffusion.m                                                                              0000640 0002626 0002066 00000002542 13122410535 014137  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [topo_new]=erosion_by_diffusion(xc,topo_old,dt_surf,dt,K)

% [topo_new]=erosion_by_diffusion(xc,topo_old,dt,K)
% J.-A. Olive, July 2015
% diffuses a topography topo_old defined on grid xc over a time dt
% using diffusivity K

N=length(xc);
topo_new=topo_old;


dt_surf=0.5*min(diff(xc)).^2/K;

if dt_surf<dt
            % DO SEVERAL TEMPERATURE SOLVES UNTIL dt is reached
            nsolve=ceil(dt/dt_surf);
            disp(['Topography erosion/deposition, nsolve=' num2str(nsolve)])
            dt_solve=dt/nsolve;
            topo_solve=topo_old;
            if (nsolve>1e5)
                disp('**Warning: topo particles too close together for erosion by diffusion');
            end;
            for ksolve=1:nsolve-1
                for i=2:N-1
                    topo_solve(i)=topo_solve(i)+dt_solve*K*(1/(0.5*(xc(i+1)+xc(i)) - 0.5*(xc(i)+xc(i-1))  ))*( (topo_solve(i+1)-topo_solve(i))/(xc(i+1)-xc(i))  -  (topo_solve(i)-topo_solve(i-1))/(xc(i)-xc(i-1))  );  
                end
                topo_solve(1)=topo_solve(2);           
            end
             topo_new=topo_solve;
else
      for i=2:N-1
        topo_new(i)=topo_old(i)+dt*K*(1/(0.5*(xc(i+1)+xc(i)) - 0.5*(xc(i)+xc(i-1))  ))*( (topo_old(i+1)-topo_old(i))/(xc(i+1)-xc(i))  -  (topo_old(i)-topo_old(i-1))/(xc(i)-xc(i-1))  );  
    end
    topo_new(1)=topo_new(2);                      
end
                                                                                                                                                              SiStER_advect_markers.m                                                                             0000640 0002626 0002066 00000002637 13122410535 014111  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [xm_new,ym_new] = SiStER_advect_markers(x,y,xm,ym,dx,dy,tstep,vx,vy)
% [xm_new,ym_new,vx_eff,vy_eff] = SiStER_advect_markers(x,y,xm,ym,dx,dy,tstep,vx,vy)
%
% Uses a Fourth-Order Runge-Kutta method to advect markers 
% in the current velocity field
%
% First cut J.-A. Olive, 2011-2012 - modified by B.Z. Klein 2013-2014


% point A (particle pts)
XA = xm;
YA = ym;                          
[qd,icn,jcn] = SiStER_locate_markers_in_grid(XA,YA,x,y,dx,dy);
[VxA,VyA] = SiStER_get_marker_velocities(qd,icn,jcn,x,y,xm,ym,vx,vy,dx,dy);

% point B
XB = XA + 0.5*tstep*VxA;
YB = YA + 0.5*tstep*VyA;
[qd,icn,jcn] = SiStER_locate_markers_in_grid(XB,YB,x,y,dx,dy);
[VxB,VyB] = SiStER_get_marker_velocities(qd,icn,jcn,x,y,xm,ym,vx,vy,dx,dy);

% point C
XC = XA + 0.5*tstep*VxB;
YC = YA + 0.5*tstep*VyB;
[qd,icn,jcn] = SiStER_locate_markers_in_grid(XC,YC,x,y,dx,dy);
[VxC,VyC] = SiStER_get_marker_velocities(qd,icn,jcn,x,y,xm,ym,vx,vy,dx,dy);

% point D
XD = XA + tstep*VxC;
YD = YA + tstep*VyC;
[qd,icn,jcn] = SiStER_locate_markers_in_grid(XD,YD,x,y,dx,dy);
[VxD,VyD] = SiStER_get_marker_velocities(qd,icn,jcn,x,y,xm,ym,vx,vy,dx,dy);


% effective velocity
vx_eff = (1/6)*(VxA + 2*VxB + 2*VxC + VxD);
vy_eff = (1/6)*(VyA + 2*VyB + 2*VyC + VyD);


%% Calculate new coordinates

PP(:,1) = XA + tstep*vx_eff;
PP(:,2) = YA + tstep*vy_eff;

xm_new = PP(:,1)';
ym_new = PP(:,2)';











                                                                                                 SiStER_assemble_L_R.m                                                                               0000640 0002626 0002066 00000051523 13275716677 013475  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [L, R, Kc, Kb]=SiStER_assemble_L_R(dx,dy,etas,etan,rho,BC,PARAMS,srhs_xx,srhs_xy,Rn,psurf)


%% Fill LHS and RHS Solution Matrices

p0cell=PARAMS.p0cell;

Nx=length(dx)+1;
Ny=length(dy)+1;

gx=PARAMS.gx;
gy=PARAMS.gy;

% Vector of right part initialization
R=zeros(Nx*Ny*3,1);
Lii = zeros(10*Nx*Ny*3,1);
Ljj = zeros(10*Nx*Ny*3,1);
Lvv = zeros(10*Nx*Ny*3,1);
nn = 1;

%%%%%%%%%% SCALING THE FD MATRIX
% Computing Kc and Kb coefficients
meta=min(min(etas));
mdx=max(dx);
mdy=max(dy);
Kc=2*meta/(mdx+mdy);
Kb=4*meta/(mdx+mdy)^2;


% pressure anchor
%IP=2;
%JP=3;
JP=ceil((Nx-1)/2); %G.Ito
IP=2;


% fill in FD matrix and right-hand side

for j=1:Nx
    for i=1:Ny
    %for j=1:Nx
        
        in=(j-1)*Ny+i;
        inp=3*in-2;
        invx=3*in-1;
        invy=3*in;
        
        
        %%%%%%% CONTINUITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        if (i==1) || (j==1) || (i==2 && j==2) || (i==Ny && j==2) || ....
           (i==2 && j==Nx) || (i==Ny && j==Nx) || (i==IP && j==JP && BC.top(2)~=3) || ...
           (BC.top(2)==3 && i==2 && j > 2 && j < Nx) %<--G.Ito
            
            
            % boundary conditions
            
            if (i==1) || (j==1) % P(i,j)=0
                %                L(inp,inp)=Kb;
                
                Lii(nn) = inp;
                Ljj(nn) = inp;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                R(inp,1)=0;
            end
            if (i==2 && j==2) || (i==Ny && j==2)  % left corners P(i,j)-P(i,j+1)=0
                %                L(inp,inp)=Kb;
                %                L(inp,inp+3*Ny)=-Kb;
                
                Lii(nn) = inp;
                Ljj(nn) = inp;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                Lii(nn) = inp;
                Ljj(nn) = inp+3*Ny;
                Lvv(nn)  = -Kb;
                nn = nn+1;
                
                R(inp,1)=0;
            end
            if (i==2 && j==Nx) || (i==Ny && j==Nx)  % right corners P(i,j)-P(i,j-1)=0
                %                L(inp,inp)=Kb;
                %                L(inp,inp-3*Ny)=-Kb;
                
                Lii(nn) = inp;
                Ljj(nn) = inp;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                Lii(nn) = inp;
                Ljj(nn) = inp-3*Ny;
                Lvv(nn)  = -Kb;
                nn = nn+1;
                
                R(inp,1)=0;
            end
            
            if (BC.top(2)==3 && i==2 && j > 2 && j < Nx)  %open top effect all nodes but corners, G.Ito
                % Pressure gradient between top two rows extrapolates to 0 pressure at very top
% 
                Lii(nn) = inp;
                Ljj(nn) = inp;
                Lvv(nn)  = Kb;
                nn = nn+1;

                R(inp,1)=Kb*psurf(j)/Kc;
               
            end
            
            if (i==IP && j==JP && BC.top(2)~=3)   % additional pressure cell i=2, j=3 P(i,j)=p0cell
                 %                L(inp,inp)=Kb;
                
                Lii(nn) = inp;
                Ljj(nn) = inp;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                R(inp,1)=Kb*p0cell/Kc;
                           
            end
           
            
        else
            
            
            % internal nodes
            
            % coeffs for vx
            %             L(inp,invx-3)=Kc/dx(j-1);
            %             L(inp,invx-3*Ny-3)=-Kc/dx(j-1);
            
            Lii(nn) = inp;
            Ljj(nn) = invx-3;
            Lvv(nn)  = Kc/dx(j-1);
            nn = nn+1;
            
            Lii(nn) = inp;
            Ljj(nn) = invx-3*Ny-3;
            Lvv(nn)  = -Kc/dx(j-1);
            nn = nn+1;
            
            % coeffs for vy
            %             L(inp,invy-3*Ny)=Kc/dy(i-1);
            %             L(inp,invy-3*Ny-3)=-Kc/dy(i-1);
            
            Lii(nn) = inp;
            Ljj(nn) = invy-3*Ny;
            Lvv(nn)  = Kc/dy(i-1);
            nn = nn+1;
            
            Lii(nn) = inp;
            Ljj(nn) = invy-3*Ny-3;
            Lvv(nn)  = -Kc/dy(i-1);
            nn = nn+1;
            
            
            % RHS
            R(inp,1)=Kc.*Rn(i,j);  %Rn for Dilation
            
        end
        
        
        %%%%%%% X-STOKES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        if ((j==1) && (i<=Ny-1)) || ((j==Nx) && (i<=Ny-1)) || (i==1 && j<=Nx-1 && j>=2) || (i==Ny-1 && j<=Nx-1 && j>=2) || (i==Ny)
            
            
            
            
            if ((j==1) && (i<=Ny-1)) %|| ((j==Nx) && (i<=Ny-1))  % X-STOKES:  left boundary 
                %              L(invx,invx)=Kb;
                
                Lii(nn) = invx;
                Ljj(nn) = invx;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                R(invx,1)=Kb*BC.left(3);
            end
            
            if ((j==Nx) && (i<=Ny-1))  %  X-STOKES:  right  boundary
                %              L(invx,invx)=Kb;
                
                Lii(nn) = invx;
                Ljj(nn) = invx;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                R(invx,1)=Kb*BC.right(3);
            end
            
            
            if i==1 && j<=Nx-1 && j>=2 % X-STOKES:  upper boundary not including sides
                if BC.top(1)==1 % free slip
                    %                  L(invx,invx+3)=Kc/(0.5*(dy(i)+dy(i+1)));
                    %                  L(invx,invx)=-Kc/(0.5*(dy(i)+dy(i+1)));
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx+3;
                    Lvv(nn)  = Kc/(0.5*(dy(i)+dy(i+1)));
                    nn = nn+1;
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx;
                    Lvv(nn)  = -Kc/(0.5*(dy(i)+dy(i+1)));
                    nn = nn+1;
                    
                    R(invx,1)=0;
                elseif BC.top(1)==0  % no slip
                    %                  L(invx,invx+3)=Kc/(dy(i)+dy(i+1));
                    %                  L(invx,invx)=-Kc*(1/dy(i)+1/(dy(i)+dy(i+1)));
                    % Multiply whole eq. by 2 so coefficients are closer to
                    % the others
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx+3;
                    Lvv(nn)  = -2*Kc/(dy(i)+dy(i+1));
                    nn = nn+1;
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx;
                    Lvv(nn)  = 2*Kc*(1/dy(i)+1/(dy(i)+dy(i+1)));
                    nn = nn+1;
                    
                    R(invx,1)=2*Kc*(1/dy(i))*BC.top_profile(j);  %G.Ito x-velocity at top boundary;
                end
                
            end
            
            if i==Ny-1 && j<=Nx-1 && j>=2 % X-STOKES: lower boundary not including sides
                if BC.bot(1)==0 % no slip
                    %                  L(invx,invx)=Kc*(1/dy(i)+1/(dy(i-1)+dy(i)));
                    %                  L(invx,invx-3)=-Kc/(dy(i-1)+dy(i));
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx;
                    Lvv(nn)  = 2*Kc*(1/dy(i)+1/(dy(i-1)+dy(i)));
                    nn = nn+1;
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx-3;
                    Lvv(nn)  = -2*Kc/(dy(i-1)+dy(i));
                    nn = nn+1;
                    
                    R(invx,1)=2*Kc*(1/dy(i))*BC.bot_profile(j);
                    
                elseif BC.bot(1)==1 % free slip
                    %                  L(invx,invx)=Kc/(0.5*(dy(i-1)+dy(i)));
                    %                  L(invx,invx-3)=-Kc/(0.5*(dy(i-1)+dy(i)));
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx;
                    Lvv(nn)  = Kc/(0.5*(dy(i-1)+dy(i)));
                    nn = nn+1;
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx-3;
                    Lvv(nn)  = -Kc/(0.5*(dy(i-1)+dy(i)));
                    nn = nn+1;
                    
                    R(invx,1)=0;
                end
                
                
                
            end
            
            if i==Ny  % lower boundary, ghost vx=0
                %              L(invx,invx)=Kb;
                
                Lii(nn) = invx;
                Ljj(nn) = invx;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                R(invx,1)=0;
            end
            
            
            
        else

            % X-STOKES: internal nodes
            
            % X-STOKES: coeffs for vx
            
            %         L(invx,invx-3*Ny)=4*etan(i+1,j)/(dx(j-1)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invx-3*Ny;
            Lvv(nn)  = 4*etan(i+1,j)/(dx(j-1)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invx,invx)=(-4/(dx(j-1)+dx(j)))*(etan(i+1,j+1)/dx(j)+etan(i+1,j)/dx(j-1))-(2/dy(i))*(etas(i+1,j)/(dy(i)+dy(i+1))+etas(i,j)/(dy(i-1)+dy(i))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invx;
            Lvv(nn)  = (-4/(dx(j-1)+dx(j)))*(etan(i+1,j+1)/dx(j)+etan(i+1,j)/dx(j-1))-(2/dy(i))*(etas(i+1,j)/(dy(i)+dy(i+1))+etas(i,j)/(dy(i-1)+dy(i)));
            nn = nn+1;
            
            %         L(invx,invx+3*Ny)=4*etan(i+1,j+1)/(dx(j)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invx+3*Ny;
            Lvv(nn)  = 4*etan(i+1,j+1)/(dx(j)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invx,invx-3)=2*etas(i,j)/(dy(i)*(dy(i-1)+dy(i))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invx-3;
            Lvv(nn)  = 2*etas(i,j)/(dy(i)*(dy(i-1)+dy(i)));
            nn = nn+1;
            
            %         L(invx,invx+3)=2*etas(i+1,j)/(dy(i)*(dy(i)+dy(i+1))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invx+3;
            Lvv(nn)  = 2*etas(i+1,j)/(dy(i)*(dy(i)+dy(i+1)));
            nn = nn+1;
            
            % X-STOKES: coeffs for vy
            
            %         L(invx,invy+3)=2*etas(i+1,j)/(dy(i)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invy+3;
            Lvv(nn)  = 2*etas(i+1,j)/(dy(i)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invx,invy+3-3*Ny)=-2*etas(i+1,j)/(dy(i)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invy+3-3*Ny;
            Lvv(nn)  = -2*etas(i+1,j)/(dy(i)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invx,invy)=-2*etas(i,j)/(dy(i)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invy;
            Lvv(nn)  = -2*etas(i,j)/(dy(i)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invx,invy-3*Ny)=2*etas(i,j)/(dy(i)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invy-3*Ny;
            Lvv(nn)  = 2*etas(i,j)/(dy(i)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            % X-STOKES: coeffs for Pp
            
            %         L(invx,inp+3*Ny+3)=-2*Kc/(dx(j-1)+dx(j)); %
            
            Lii(nn) = invx;
            Ljj(nn) = inp+3*Ny+3;
            Lvv(nn)  = -2*Kc/(dx(j-1)+dx(j));
            nn = nn+1;
            
            %         L(invx,inp+3)=2*Kc/(dx(j-1)+dx(j)); %
            
            Lii(nn) = invx;
            Ljj(nn) = inp+3;
            Lvv(nn)  = 2*Kc/(dx(j-1)+dx(j));
            nn = nn+1;
            
            % RHS
            R(invx,1)=-gx*(rho(i,j)+rho(i+1,j))/2 - ...
                2.*(srhs_xx(i+1,j+1)-srhs_xx(i+1,j))/(dx(j)+dx(j-1)) - ...
                   (srhs_xy(i+1,j)-srhs_xy(i,j))/dy(i) + ...
                2*(etan(i+1,j+1).*Rn(i+1,j+1) -...
                    etan(i+1,j).*Rn(i+1,j))./(dx(j)+dx(j-1));  %dilation with plasticity
            
        end
        
        
        
        
        %%%%%%% Y-STOKES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        if (j==1 && i<=Ny-1 && i>=2) || (j==Nx-1 && i<=Ny-1 && i>=2) || (j==Nx) || (i==1 && j<=Nx-1) || (i==Ny && j<=Nx-1)
                   
            if j==1 && i>=2 && i<=Ny-1  % left boundary without edges
                if BC.left(1)==1 % free slip
                    %              L(invy,invy)=-2*Kc/(dx(j)+dx(j+1));
                    %              L(invy,invy+3*Ny)=2*Kc/(dx(j)+dx(j+1));
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn)  = -2*Kc/(dx(j)+dx(j+1));
                    nn = nn+1;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy+3*Ny;
                    Lvv(nn)  = 2*Kc/(dx(j)+dx(j+1));
                    nn = nn+1;
                    
                    R(invy,1)=0;
                    
                elseif BC.left(1)==0 % no slip
                    %              L(invy,invy)=-2*Kc*(1/(dx(j)+dx(j+1))+1/dx(j));
                    %              L(invy,invy+3*Ny)=2*Kc/(dx(j+1)+dx(j));
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn)  = -2*Kc*(1/(dx(j)+dx(j+1))+1/dx(j));
                    nn = nn+1;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy+3*Ny;
                    Lvv(nn)  = 2*Kc/(dx(j+1)+dx(j));
                    nn = nn+1;
                    
                    R(invy,1)=0;
                    
                end
            end
            if j==Nx-1 && i<=Ny-1 && i>=2  % right boundary without edges 
                if BC.right(1)==1 % free slip
                    %              L(invy,invy)=Kc/(0.5*(dx(j)+dx(j-1)));
                    %              L(invy,invy-3*Ny)=-Kc/(0.5*(dx(j)+dx(j-1)));
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn)  = Kc/(0.5*(dx(j)+dx(j-1)));
                    nn = nn+1;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy-3*Ny;
                    Lvv(nn)  = -Kc/(0.5*(dx(j)+dx(j-1)));
                    nn = nn+1;
                    
                    R(invy,1)=0;
                    
                elseif BC.right(1)==0 % no slip
                    %              L(invy,invy)=Kc*(1/(dx(j)+dx(j-1))+1/dx(j));
                    %              L(invy,invy-3*Ny)=-Kc/(dx(j)+dx(j-1));
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn)  = Kc*(1/(dx(j)+dx(j-1))+1/dx(j));
                    nn = nn+1;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy-3*Ny;
                    Lvv(nn) = -Kc/(dx(j)+dx(j-1));
                    nn = nn+1;
                    
                    R(invy,1)=0;
                    
                end
                
            end
            if j==Nx  % right boundary, ghost vy=0
                %              L(invy,invy)=Kb;
                
                Lii(nn) = invy;
                Ljj(nn) = invy;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                R(invy,1)=0;
                
            end
            
            if i==1 && j<=Nx-1  % upper boundary
                if BC.top(2)==0
                    %              L(invy,invy)=Kb;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn)  = Kb;
                    nn = nn+1;
                    
                    R(invy,1)=Kb*BC.top(3);
                elseif BC.top(2)==3
                    %               open top dvy/dy=0, G.Ito
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn) = Kb;
                    nn      = nn+1;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy+3;
                    Lvv(nn) = -Kb;
                    nn      = nn+1;
                    
                    R(invy,1)=0;
                    
                end
            end
            
            if i==Ny && j<=Nx-1  % lower boundary without edges
                
                if BC.bot(2)==0
                    %              L(invy,invy)=Kb;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn)  = Kb;
                    nn = nn+1;
                    
                    R(invy,1)=Kb*BC.bot(3);
                    
                elseif BC.bot(2) == 2 %%  external boundary
                    
                    dL=300e3;
                                        
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn) = Kc*(1/dL+1/dy(i-1));
                    nn = nn+1;
                    
                    Lii(nn) = invy; 
                    Ljj(nn) = invy-3;
                    Lvv(nn) = Kc*(-1/dy(i-1));
                    nn = nn+1;

                    R(invy, 1) = 0;
                end
                    
            end
            
            
        else
            
            
            
            % internal nodes
            
            % coeffs for vy
            %         L(invy,invy-3*Ny)=2*etas(i,j)/(dx(j)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invy-3*Ny;
            Lvv(nn)  = 2*etas(i,j)/(dx(j)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invy,invy)=(-4/(dy(i-1)+dy(i)))*(etan(i+1,j+1)/dy(i)+etan(i,j+1)/dy(i-1))-(2/dx(j))*(etas(i,j+1)/(dx(j)+dx(j+1))+etas(i,j)/(dx(j-1)+dx(j))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invy;
            Lvv(nn)  = (-4/(dy(i-1)+dy(i)))*(etan(i+1,j+1)/dy(i)+etan(i,j+1)/dy(i-1))-(2/dx(j))*(etas(i,j+1)/(dx(j)+dx(j+1))+etas(i,j)/(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invy,invy+3*Ny)=2*etas(i,j+1)/(dx(j)*(dx(j)+dx(j+1))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invy+3*Ny;
            Lvv(nn)  = 2*etas(i,j+1)/(dx(j)*(dx(j)+dx(j+1)));
            nn = nn+1;
            
            %         L(invy,invy-3)=4*etan(i,j+1)/(dy(i-1)*(dy(i-1)+dy(i))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invy-3;
            Lvv(nn)  = 4*etan(i,j+1)/(dy(i-1)*(dy(i-1)+dy(i)));
            nn  = nn+1;
            
            %         L(invy,invy+3)=4*etan(i+1,j+1)/(dy(i)*(dy(i-1)+dy(i))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invy+3;
            Lvv(nn)  = 4*etan(i+1,j+1)/(dy(i)*(dy(i-1)+dy(i)));
            nn = nn+1;
            
            % coeffs for vx
            
            %         L(invy,invx+3*Ny)=2*etas(i,j+1)/(dx(j)*(dy(i-1)+dy(i)));  %
            
            Lii(nn) = invy;
            Ljj(nn) = invx+3*Ny;
            Lvv(nn)  = 2*etas(i,j+1)/(dx(j)*(dy(i-1)+dy(i)));
            nn = nn+1;
            
            %         L(invy,invx+3*Ny-3)=-2*etas(i,j+1)/(dx(j)*(dy(i-1)+dy(i))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invx+3*Ny-3;
            Lvv(nn)  = -2*etas(i,j+1)/(dx(j)*(dy(i-1)+dy(i)));
            nn = nn+1;
            
            
            %         L(invy,invx)=-2*etas(i,j)/(dx(j)*(dy(i-1)+dy(i))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invx;
            Lvv(nn)  = -2*etas(i,j)/(dx(j)*(dy(i-1)+dy(i)));
            nn = nn+1;
            
            %         L(invy,invx-3)=2*etas(i,j)/(dx(j)*(dy(i-1)+dy(i))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invx-3;
            Lvv(nn)  = 2*etas(i,j)/(dx(j)*(dy(i-1)+dy(i)));
            nn=nn+1;
            
            % coeffs for Pp
            %         L(invy,inp+3*Ny+3)=-2*Kc/(dy(i-1)+dy(i)); %
            
            Lii(nn) = invy;
            Ljj(nn) = inp+3*Ny+3;
            Lvv(nn)  = -2*Kc/(dy(i-1)+dy(i));
            nn=nn+1;
            
            %         L(invy,inp+3*Ny)=2*Kc/(dy(i-1)+dy(i)); %
            
            Lii(nn) = invy;
            Ljj(nn) = inp+3*Ny;
            Lvv(nn)  = 2*Kc/(dy(i-1)+dy(i));
            nn=nn+1;
            
            % RHS (note: sxx = -syy because deviatoric, but not with
            % dilation!)
            R(invy,1)=-gy*(rho(i,j)+rho(i,j+1))/2 + ...
              2.*(srhs_xx(i+1,j+1)-srhs_xx(i,j+1))/(dy(i)+dy(i-1)) - ...
                 (srhs_xy(i,j+1)-srhs_xy(i,j))/dx(j) +...
              2*(etan(i+1,j+1).*Rn(i+1,j+1) -...
                  etan(i,j+1).*Rn(i,j+1))./(dy(i)+dy(i-1));%
            
        end
        

        
    end % for j
end % for i

nn = nn-1;
%%% Build Sparse Matrix

Li = Lii(1:nn);
Lj = Ljj(1:nn);
Lv = Lvv(1:nn);


L = sparse(Li,Lj,Lv);
% SOLVE LINEAR SYSTEM now down in SiStER_flow_solve G.Ito
% Matlab direct solver
% tic
% S=L\R;
% toc
% 





end

                                                                                                                                                                             SiStER_assemble_L_R_new.m                                                                           0000640 0002626 0002066 00000053150 13225074633 014324  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [L, R, Kc, Kb]=SiStER_assemble_L_R_Olive(dx,dy,etas,etan,rho,BC,PARAMS,srhs_xx,srhs_xy,psurf)

%% BOUNDARY CONDITIONS 

% Flow boundary conditions
% BC....(1): tangential velocity
%    0 = No slip
%    1 = Free slip (if using multiple side BCs, set in BC.side, not
%        BC.side_w/l/m
% BC...(2): velocity normal to boundary 
%    0 = constant, conservative 
%    1 = triangular profile (top/bot only)
%    2 = multiple (sides only)
%    3 = open     (only for top, bottom, and multiple)
% BC...(3): value of normal velocity (or max. value along the profile, if BC..(2) is 1)
%
% BC.bot(4): value of tangential velocity for No Slip (only works for bottom boundary)

%% Fill LHS and RHS Solution Matrices

p0cell=PARAMS.p0cell;

Nx=length(dx)+1;
Ny=length(dy)+1;

gx=PARAMS.gx;
gy=PARAMS.gy;

% Vector of right part initialization
R=zeros(Nx*Ny*3,1);
Lii = zeros(10*Nx*Ny*3,1);
Ljj = zeros(10*Nx*Ny*3,1);
Lvv = zeros(10*Nx*Ny*3,1);
nn = 1;

%%%%%%%%%% SCALING THE FD MATRIX
% Computing Kc and Kb coefficients
meta=min(min(etas));
mdx=max(dx);
mdy=max(dy);
Kc=2*meta/(mdx+mdy);
Kb=4*meta/(mdx+mdy)^2;


% pressure anchor
%IP=2;
%JP=3;
JP=ceil((Nx-1)/2); %G.Ito
JP=Nx; %G.Ito debugging
IP=2;


% fill in FD matrix and right-hand side

for j=1:Nx
    for i=1:Ny
    %for j=1:Nx
        
        in=(j-1)*Ny+i;
        inp=3*in-2;
        invx=3*in-1;
        invy=3*in;
        
        
        %%%%%%% CONTINUITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
%         if (i==1) || (j==1) || (i==2 && j==2) || (i==2 && j==Nx) || ....
%            (i==Ny && j==2) || (i==Ny && j==Nx) || (i==IP && j==JP && BC.top(2)~=3) || ...
%            (BC.top(2)==3 && i==2 && j > 2 && j < Ny) %<--G.Ito

        if (i==1) || (j==1) || (i==2 && j==2) || (i==Ny && j==2) || ...
           (i==Ny && j==Nx) || (i==IP && j==JP) || ...
           (BC.top(2)==3 && i==2 && j > 2 && j < Nx) %debugging
            
            
            % boundary conditions
            
            if (i==1) || (j==1) % P(i,j)=0
                %                L(inp,inp)=Kb;
                
                Lii(nn) = inp;
                Ljj(nn) = inp;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                R(inp,1)=0;
            end
            if (i==2 && j==2) || (i==Ny && j==2)  % left corners P(i,j)-P(i,j+1)=0
                %                L(inp,inp)=Kb;
                %                L(inp,inp+3*Ny)=-Kb;
                
                Lii(nn) = inp;
                Ljj(nn) = inp;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                Lii(nn) = inp;
                Ljj(nn) = inp+3*Ny;
                Lvv(nn)  = -Kb;
                nn = nn+1;
                
                R(inp,1)=0;
            end
            %if (i==2 && j==Nx) || (i==Ny && j==Nx)  % right corners P(i,j)-P(i,j-1)=0
            if (i==Ny && j==Nx)  % debugging
               %                L(inp,inp)=Kb;
                %                L(inp,inp-3*Ny)=-Kb;
                
                Lii(nn) = inp;
                Ljj(nn) = inp;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                Lii(nn) = inp;
                Ljj(nn) = inp-3*Ny;
                Lvv(nn)  = -Kb;
                nn = nn+1;
                
                R(inp,1)=0;
            end
            
            %if (BC.top(2)==3 && i==2 && j > 2 && j < Ny)  %open top effects all nodes but corners, G.Ito
                           % Pressure gradient between top two rows extrapolates to 0 pressure at very top
            if (BC.top(2)==3 && i==2 && j > 2 && j < Nx)  %debugging
 % 
                Lii(nn) = inp;
                Ljj(nn) = inp;
                Lvv(nn)  = Kb.*(1+dy(i-1)/(dy(i)+dy(i-1)));
                nn = nn+1;

                Lii(nn) = inp;
                Ljj(nn) = inp+3;
                Lvv(nn)  = -Kb.*dy(i-1)/(dy(i)+dy(i-1));
                nn = nn+1;

                R(inp,1)=psurf(j);
               
            end
            
            %if (i==IP && j==JP && BC.top(2)~=3)   % additional pressure cell i=2, j=3 P(i,j)=p0cell
            if (i==IP && j==JP)   % debugging
                 %                L(inp,inp)=Kb;
                
                Lii(nn) = inp;
                Ljj(nn) = inp;
                Lvv(nn)  = Kb;
                nn = nn+1;
                R(inp,1)=Kb*p0cell/Kc;
            end
           
            
        else
            
            
            % internal nodes
            
            % coeffs for vx
            %             L(inp,invx-3)=Kc/dx(j-1);
            %             L(inp,invx-3*Ny-3)=-Kc/dx(j-1);
            
            Lii(nn) = inp;
            Ljj(nn) = invx-3;
            Lvv(nn)  = Kc/dx(j-1);
            nn = nn+1;
            
            Lii(nn) = inp;
            Ljj(nn) = invx-3*Ny-3;
            Lvv(nn)  = -Kc/dx(j-1);
            nn = nn+1;
            
            % coeffs for vy
            %             L(inp,invy-3*Ny)=Kc/dy(i-1);
            %             L(inp,invy-3*Ny-3)=-Kc/dy(i-1);
            
            Lii(nn) = inp;
            Ljj(nn) = invy-3*Ny;
            Lvv(nn)  = Kc/dy(i-1);
            nn = nn+1;
            
            Lii(nn) = inp;
            Ljj(nn) = invy-3*Ny-3;
            Lvv(nn)  = -Kc/dy(i-1);
            nn = nn+1;
            
            
            % RHS
            R(inp,1)=0;%-drhodt(i,j);
            
        end
        
        
        %%%%%%% X-STOKES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        if ((j==1) && (i<=Ny-1)) || ((j==Nx) && (i<=Ny-1)) || (i==1 && j<=Nx-1 && j>=2) || (i==Ny-1 && j<=Nx-1 && j>=2) || (i==Ny)
            
            
            
            
            if ((j==1) && (i<=Ny-1)) %|| ((j==Nx) && (i<=Ny-1))  % X-STOKES:  left boundary 
                %              L(invx,invx)=Kb;
                
                Lii(nn) = invx;
                Ljj(nn) = invx;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                R(invx,1)=Kb*BC.left(3);
            end
            
            if ((j==Nx) && (i<=Ny-1))  %  X-STOKES:  right  boundary
                %              L(invx,invx)=Kb;
                
                Lii(nn) = invx;
                Ljj(nn) = invx;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                R(invx,1)=Kb*BC.right(3);
            end
            
            
            if i==1 && j<=Nx-1 && j>=2 % X-STOKES:  upper boundary not including sides
                if BC.top(1)==1 % free slip
                    %                  L(invx,invx+3)=Kc/(0.5*(dy(i)+dy(i+1)));
                    %                  L(invx,invx)=-Kc/(0.5*(dy(i)+dy(i+1)));
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx+3;
                    Lvv(nn)  = Kc/(0.5*(dy(i)+dy(i+1)));
                    nn = nn+1;
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx;
                    Lvv(nn)  = -Kc/(0.5*(dy(i)+dy(i+1)));
                    nn = nn+1;
                    
                    R(invx,1)=0;
                elseif BC.top(1)==0  % no slip
                    %                  L(invx,invx+3)=Kc/(dy(i)+dy(i+1));
                    %                  L(invx,invx)=-Kc*(1/dy(i)+1/(dy(i)+dy(i+1)));
                    % Multiply whole eq. by 2 so coefficients are closer to
                    % the others
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx+3;
                    Lvv(nn)  = -2*Kc/(dy(i)+dy(i+1));
                    nn = nn+1;
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx;
                    Lvv(nn)  = 2*Kc*(1/dy(i)+1/(dy(i)+dy(i+1)));
                    nn = nn+1;
                    
                    R(invx,1)=2*Kc*(1/dy(i))*BC.top_profile(j);  %G.Ito x-velocity at top boundary;
                end
                
            end
            
            if i==Ny-1 && j<=Nx-1 && j>=2 % X-STOKES: lower boundary not including sides
                if BC.bot(1)==0 % no slip
                    %                  L(invx,invx)=Kc*(1/dy(i)+1/(dy(i-1)+dy(i)));
                    %                  L(invx,invx-3)=-Kc/(dy(i-1)+dy(i));
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx;
                    Lvv(nn)  = 2*Kc*(1/dy(i)+1/(dy(i-1)+dy(i)));
                    nn = nn+1;
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx-3;
                    Lvv(nn)  = -2*Kc/(dy(i-1)+dy(i));
                    nn = nn+1;
                    
                    R(invx,1)=2*Kc*(1/dy(i))*BC.bot_profile(j);
                    
                elseif BC.bot(1)==1 % free slip
                    %                  L(invx,invx)=Kc/(0.5*(dy(i-1)+dy(i)));
                    %                  L(invx,invx-3)=-Kc/(0.5*(dy(i-1)+dy(i)));
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx;
                    Lvv(nn)  = Kc/(0.5*(dy(i-1)+dy(i)));
                    nn = nn+1;
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx-3;
                    Lvv(nn)  = -Kc/(0.5*(dy(i-1)+dy(i)));
                    nn = nn+1;
                    
                    R(invx,1)=0;
                end
                
                
                
            end
            
            if i==Ny  % lower boundary, ghost vx=0
                %              L(invx,invx)=Kb;
                
                Lii(nn) = invx;
                Ljj(nn) = invx;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                R(invx,1)=0;
            end
            
            
            
        else

            % X-STOKES: internal nodes
            
            % X-STOKES: coeffs for vx
            
            %         L(invx,invx-3*Ny)=4*etan(i+1,j)/(dx(j-1)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invx-3*Ny;
            Lvv(nn)  = 4*etan(i+1,j)/(dx(j-1)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invx,invx)=(-4/(dx(j-1)+dx(j)))*(etan(i+1,j+1)/dx(j)+etan(i+1,j)/dx(j-1))-(2/dy(i))*(etas(i+1,j)/(dy(i)+dy(i+1))+etas(i,j)/(dy(i-1)+dy(i))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invx;
            Lvv(nn)  = (-4/(dx(j-1)+dx(j)))*(etan(i+1,j+1)/dx(j)+etan(i+1,j)/dx(j-1))-(2/dy(i))*(etas(i+1,j)/(dy(i)+dy(i+1))+etas(i,j)/(dy(i-1)+dy(i)));
            nn = nn+1;
            
            %         L(invx,invx+3*Ny)=4*etan(i+1,j+1)/(dx(j)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invx+3*Ny;
            Lvv(nn)  = 4*etan(i+1,j+1)/(dx(j)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invx,invx-3)=2*etas(i,j)/(dy(i)*(dy(i-1)+dy(i))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invx-3;
            Lvv(nn)  = 2*etas(i,j)/(dy(i)*(dy(i-1)+dy(i)));
            nn = nn+1;
            
            %         L(invx,invx+3)=2*etas(i+1,j)/(dy(i)*(dy(i)+dy(i+1))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invx+3;
            Lvv(nn)  = 2*etas(i+1,j)/(dy(i)*(dy(i)+dy(i+1)));
            nn = nn+1;
            
            % X-STOKES: coeffs for vy
            
            %         L(invx,invy+3)=2*etas(i+1,j)/(dy(i)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invy+3;
            Lvv(nn)  = 2*etas(i+1,j)/(dy(i)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invx,invy+3-3*Ny)=-2*etas(i+1,j)/(dy(i)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invy+3-3*Ny;
            Lvv(nn)  = -2*etas(i+1,j)/(dy(i)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invx,invy)=-2*etas(i,j)/(dy(i)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invy;
            Lvv(nn)  = -2*etas(i,j)/(dy(i)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invx,invy-3*Ny)=2*etas(i,j)/(dy(i)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invy-3*Ny;
            Lvv(nn)  = 2*etas(i,j)/(dy(i)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            % X-STOKES: coeffs for Pp
            
            %         L(invx,inp+3*Ny+3)=-2*Kc/(dx(j-1)+dx(j)); %
            
            Lii(nn) = invx;
            Ljj(nn) = inp+3*Ny+3;
            Lvv(nn)  = -2*Kc/(dx(j-1)+dx(j));
            nn = nn+1;
            
            %         L(invx,inp+3)=2*Kc/(dx(j-1)+dx(j)); %
            
            Lii(nn) = invx;
            Ljj(nn) = inp+3;
            Lvv(nn)  = 2*Kc/(dx(j-1)+dx(j));
            nn = nn+1;
            
            % RHS
            R(invx,1)=-gx*(rho(i,j)+rho(i+1,j))/2 - ...
                2.*(srhs_xx(i+1,j+1)-srhs_xx(i+1,j))/(dx(j)+dx(j-1)) - ...
                (srhs_xy(i+1,j)-srhs_xy(i,j))/dy(i);
            
        end
        
        
        
        
        %%%%%%% Y-STOKES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        if (j==1 && i<=Ny-1 && i>=2) || (j==Nx-1 && i<=Ny-1 && i>=2) || (j==Nx) || (i==1 && j<=Nx-1) || (i==Ny && j<=Nx-1)
                   
            if j==1 && i>=2 && i<=Ny-1  % left boundary without edges
                if BC.left(1)==1 % free slip
                    %              L(invy,invy)=-2*Kc/(dx(j)+dx(j+1));
                    %              L(invy,invy+3*Ny)=2*Kc/(dx(j)+dx(j+1));
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn)  = -2*Kc/(dx(j)+dx(j+1));
                    nn = nn+1;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy+3*Ny;
                    Lvv(nn)  = 2*Kc/(dx(j)+dx(j+1));
                    nn = nn+1;
                    
                    R(invy,1)=0;
                    
                elseif BC.left(1)==0 % no slip
                    %              L(invy,invy)=-2*Kc*(1/(dx(j)+dx(j+1))+1/dx(j));
                    %              L(invy,invy+3*Ny)=2*Kc/(dx(j+1)+dx(j));
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn)  = -2*Kc*(1/(dx(j)+dx(j+1))+1/dx(j));
                    nn = nn+1;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy+3*Ny;
                    Lvv(nn)  = 2*Kc/(dx(j+1)+dx(j));
                    nn = nn+1;
                    
                    R(invy,1)=0;
                    
                end
            end
            if j==Nx-1 && i<=Ny-1 && i>=2  % right boundary without edges 
                if BC.right(1)==1 % free slip
                    %              L(invy,invy)=Kc/(0.5*(dx(j)+dx(j-1)));
                    %              L(invy,invy-3*Ny)=-Kc/(0.5*(dx(j)+dx(j-1)));
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn)  = Kc/(0.5*(dx(j)+dx(j-1)));
                    nn = nn+1;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy-3*Ny;
                    Lvv(nn)  = -Kc/(0.5*(dx(j)+dx(j-1)));
                    nn = nn+1;
                    
                    R(invy,1)=0;
                    
                elseif BC.right(1)==0 % no slip
                    %              L(invy,invy)=Kc*(1/(dx(j)+dx(j-1))+1/dx(j));
                    %              L(invy,invy-3*Ny)=-Kc/(dx(j)+dx(j-1));
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn)  = Kc*(1/(dx(j)+dx(j-1))+1/dx(j));
                    nn = nn+1;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy-3*Ny;
                    Lvv(nn) = -Kc/(dx(j)+dx(j-1));
                    nn = nn+1;
                    
                    R(invy,1)=0;
                    
                end
                
            end
            if j==Nx  % right boundary, ghost vy=0
                %              L(invy,invy)=Kb;
                
                Lii(nn) = invy;
                Ljj(nn) = invy;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                R(invy,1)=0;
                
            end
            
            if i==1 && j<=Nx-1  % upper boundary
                if BC.top(2)==0
                    %              L(invy,invy)=Kb;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn)  = Kb;
                    nn = nn+1;
                    
                    R(invy,1)=Kb*BC.top(3);
                elseif BC.top(2)==3
                    %               open top dvy/dy=0, G.Ito
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn) = Kb;
                    nn      = nn+1;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy+3;
                    Lvv(nn) = -Kb;
                    nn      = nn+1;
                    
                    R(invy,1)=0;
                    
                end
            end
            
            if i==Ny && j<=Nx-1  % lower boundary without edges
                
                if BC.bot(2)==0
                    %              L(invy,invy)=Kb;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn)  = Kb;
                    nn = nn+1;
                    
                    R(invy,1)=Kb*BC.bot(3);
                    
                elseif BC.bot(2) == 2 %%  external boundary
                    
                    dL=300e3;
                                        
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn) = Kc*(1/dL+1/dy(i-1));
                    nn = nn+1;
                    
                    Lii(nn) = invy; 
                    Ljj(nn) = invy-3;
                    Lvv(nn) = Kc*(-1/dy(i-1));
                    nn = nn+1;

                    R(invy, 1) = 0;
                end
                    
            end
            
            
        else
            
            
            
            % internal nodes
            
            % coeffs for vy
            %         L(invy,invy-3*Ny)=2*etas(i,j)/(dx(j)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invy-3*Ny;
            Lvv(nn)  = 2*etas(i,j)/(dx(j)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invy,invy)=(-4/(dy(i-1)+dy(i)))*(etan(i+1,j+1)/dy(i)+etan(i,j+1)/dy(i-1))-(2/dx(j))*(etas(i,j+1)/(dx(j)+dx(j+1))+etas(i,j)/(dx(j-1)+dx(j))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invy;
            Lvv(nn)  = (-4/(dy(i-1)+dy(i)))*(etan(i+1,j+1)/dy(i)+etan(i,j+1)/dy(i-1))-(2/dx(j))*(etas(i,j+1)/(dx(j)+dx(j+1))+etas(i,j)/(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invy,invy+3*Ny)=2*etas(i,j+1)/(dx(j)*(dx(j)+dx(j+1))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invy+3*Ny;
            Lvv(nn)  = 2*etas(i,j+1)/(dx(j)*(dx(j)+dx(j+1)));
            nn = nn+1;
            
            %         L(invy,invy-3)=4*etan(i,j+1)/(dy(i-1)*(dy(i-1)+dy(i))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invy-3;
            Lvv(nn)  = 4*etan(i,j+1)/(dy(i-1)*(dy(i-1)+dy(i)));
            nn  = nn+1;
            
            %         L(invy,invy+3)=4*etan(i+1,j+1)/(dy(i)*(dy(i-1)+dy(i))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invy+3;
            Lvv(nn)  = 4*etan(i+1,j+1)/(dy(i)*(dy(i-1)+dy(i)));
            nn = nn+1;
            
            % coeffs for vx
            
            %         L(invy,invx+3*Ny)=2*etas(i,j+1)/(dx(j)*(dy(i-1)+dy(i)));  %
            
            Lii(nn) = invy;
            Ljj(nn) = invx+3*Ny;
            Lvv(nn)  = 2*etas(i,j+1)/(dx(j)*(dy(i-1)+dy(i)));
            nn = nn+1;
            
            %         L(invy,invx+3*Ny-3)=-2*etas(i,j+1)/(dx(j)*(dy(i-1)+dy(i))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invx+3*Ny-3;
            Lvv(nn)  = -2*etas(i,j+1)/(dx(j)*(dy(i-1)+dy(i)));
            nn = nn+1;
            
            
            %         L(invy,invx)=-2*etas(i,j)/(dx(j)*(dy(i-1)+dy(i))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invx;
            Lvv(nn)  = -2*etas(i,j)/(dx(j)*(dy(i-1)+dy(i)));
            nn = nn+1;
            
            %         L(invy,invx-3)=2*etas(i,j)/(dx(j)*(dy(i-1)+dy(i))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invx-3;
            Lvv(nn)  = 2*etas(i,j)/(dx(j)*(dy(i-1)+dy(i)));
            nn=nn+1;
            
            % coeffs for Pp
            %         L(invy,inp+3*Ny+3)=-2*Kc/(dy(i-1)+dy(i)); %
            
            Lii(nn) = invy;
            Ljj(nn) = inp+3*Ny+3;
            Lvv(nn)  = -2*Kc/(dy(i-1)+dy(i));
            nn=nn+1;
            
            %         L(invy,inp+3*Ny)=2*Kc/(dy(i-1)+dy(i)); %
            
            Lii(nn) = invy;
            Ljj(nn) = inp+3*Ny;
            Lvv(nn)  = 2*Kc/(dy(i-1)+dy(i));
            nn=nn+1;
            
            % RHS (note: sxx = -syy because deviatoric)
            R(invy,1)=-gy*(rho(i,j)+rho(i,j+1))/2 + ...
                2.*(srhs_xx(i+1,j+1)-srhs_xx(i,j+1))/(dy(i)+dy(i-1)) - ...
                (srhs_xy(i,j+1)-srhs_xy(i,j))/dx(j); %
            
        end
        

        
    end % for j
end % for i

nn = nn-1;
%%% Build Sparse Matrix

Li = Lii(1:nn);
Lj = Ljj(1:nn);
Lv = Lvv(1:nn);


L = sparse(Li,Lj,Lv);
% SOLVE LINEAR SYSTEM now down in SiStER_flow_solve G.Ito
% Matlab direct solver
% tic
% S=L\R;
% toc
% 





end

                                                                                                                                                                                                                                                                                                                                                                                                                        SiStER_assemble_L_R_orig.m                                                                          0000640 0002626 0002066 00000051304 13225075365 014475  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [L, R, Kc, Kb]=SiStER_assemble_L_R_Olive(dx,dy,etas,etan,rho,BC,PARAMS,srhs_xx,srhs_xy,psurf)


%% Fill LHS and RHS Solution Matrices

p0cell=PARAMS.p0cell;

Nx=length(dx)+1;
Ny=length(dy)+1;

gx=PARAMS.gx;
gy=PARAMS.gy;

% Vector of right part initialization
R=zeros(Nx*Ny*3,1);
Lii = zeros(10*Nx*Ny*3,1);
Ljj = zeros(10*Nx*Ny*3,1);
Lvv = zeros(10*Nx*Ny*3,1);
nn = 1;

%%%%%%%%%% SCALING THE FD MATRIX
% Computing Kc and Kb coefficients
meta=min(min(etas));
mdx=max(dx);
mdy=max(dy);
Kc=2*meta/(mdx+mdy);
Kb=4*meta/(mdx+mdy)^2;


% pressure anchor
%IP=2;
%JP=3;
JP=ceil((Nx-1)/2); %G.Ito
IP=2;


% fill in FD matrix and right-hand side

for j=1:Nx
    for i=1:Ny
    %for j=1:Nx
        
        in=(j-1)*Ny+i;
        inp=3*in-2;
        invx=3*in-1;
        invy=3*in;
        
        
        %%%%%%% CONTINUITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        if (i==1) || (j==1) || (i==2 && j==2) || (i==2 && j==Nx) || ....
           (i==Ny && j==2) || (i==Ny && j==Nx) || (i==IP && j==JP && BC.top(2)~=3) || ...
           (BC.top(2)==3 && i==2 && j > 2 && j < Ny) %<--G.Ito
            
            
            % boundary conditions
            
            if (i==1) || (j==1) % P(i,j)=0
                %                L(inp,inp)=Kb;
                
                Lii(nn) = inp;
                Ljj(nn) = inp;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                R(inp,1)=0;
            end
            if (i==2 && j==2) || (i==Ny && j==2)  % left corners P(i,j)-P(i,j+1)=0
                %                L(inp,inp)=Kb;
                %                L(inp,inp+3*Ny)=-Kb;
                
                Lii(nn) = inp;
                Ljj(nn) = inp;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                Lii(nn) = inp;
                Ljj(nn) = inp+3*Ny;
                Lvv(nn)  = -Kb;
                nn = nn+1;
                
                R(inp,1)=0;
            end
            if (i==2 && j==Nx) || (i==Ny && j==Nx)  % right corners P(i,j)-P(i,j-1)=0
                %                L(inp,inp)=Kb;
                %                L(inp,inp-3*Ny)=-Kb;
                
                Lii(nn) = inp;
                Ljj(nn) = inp;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                Lii(nn) = inp;
                Ljj(nn) = inp-3*Ny;
                Lvv(nn)  = -Kb;
                nn = nn+1;
                
                R(inp,1)=0;
            end
            
            if (BC.top(2)==3 && i==2 && j > 2 && j < Ny)  %open top effect all nodes but corners, G.Ito
                % Pressure gradient between top two rows extrapolates to 0 pressure at very top
% 
                Lii(nn) = inp;
                Ljj(nn) = inp;
                Lvv(nn)  = Kb.*(1+dy(i-1)/(dy(i)+dy(i-1)));
                nn = nn+1;

                Lii(nn) = inp;
                Ljj(nn) = inp+3;
                Lvv(nn)  = -Kb.*dy(i-1)/(dy(i)+dy(i-1));
                nn = nn+1;

                R(inp,1)=0;
               
            end
            
            if (i==IP && j==JP && BC.top(2)~=3)   % additional pressure cell i=2, j=3 P(i,j)=p0cell
                 %                L(inp,inp)=Kb;
                
                Lii(nn) = inp;
                Ljj(nn) = inp;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                R(inp,1)=Kb*p0cell/Kc;
                           
            end
           
            
        else
            
            
            % internal nodes
            
            % coeffs for vx
            %             L(inp,invx-3)=Kc/dx(j-1);
            %             L(inp,invx-3*Ny-3)=-Kc/dx(j-1);
            
            Lii(nn) = inp;
            Ljj(nn) = invx-3;
            Lvv(nn)  = Kc/dx(j-1);
            nn = nn+1;
            
            Lii(nn) = inp;
            Ljj(nn) = invx-3*Ny-3;
            Lvv(nn)  = -Kc/dx(j-1);
            nn = nn+1;
            
            % coeffs for vy
            %             L(inp,invy-3*Ny)=Kc/dy(i-1);
            %             L(inp,invy-3*Ny-3)=-Kc/dy(i-1);
            
            Lii(nn) = inp;
            Ljj(nn) = invy-3*Ny;
            Lvv(nn)  = Kc/dy(i-1);
            nn = nn+1;
            
            Lii(nn) = inp;
            Ljj(nn) = invy-3*Ny-3;
            Lvv(nn)  = -Kc/dy(i-1);
            nn = nn+1;
            
            
            % RHS
            R(inp,1)=0;%-drhodt(i,j);
            
        end
        
        
        %%%%%%% X-STOKES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        if ((j==1) && (i<=Ny-1)) || ((j==Nx) && (i<=Ny-1)) || (i==1 && j<=Nx-1 && j>=2) || (i==Ny-1 && j<=Nx-1 && j>=2) || (i==Ny)
            
            
            
            
            if ((j==1) && (i<=Ny-1)) %|| ((j==Nx) && (i<=Ny-1))  % X-STOKES:  left boundary 
                %              L(invx,invx)=Kb;
                
                Lii(nn) = invx;
                Ljj(nn) = invx;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                R(invx,1)=Kb*BC.left(3);
            end
            
            if ((j==Nx) && (i<=Ny-1))  %  X-STOKES:  right  boundary
                %              L(invx,invx)=Kb;
                
                Lii(nn) = invx;
                Ljj(nn) = invx;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                R(invx,1)=Kb*BC.right(3);
            end
            
            
            if i==1 && j<=Nx-1 && j>=2 % X-STOKES:  upper boundary not including sides
                if BC.top(1)==1 % free slip
                    %                  L(invx,invx+3)=Kc/(0.5*(dy(i)+dy(i+1)));
                    %                  L(invx,invx)=-Kc/(0.5*(dy(i)+dy(i+1)));
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx+3;
                    Lvv(nn)  = Kc/(0.5*(dy(i)+dy(i+1)));
                    nn = nn+1;
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx;
                    Lvv(nn)  = -Kc/(0.5*(dy(i)+dy(i+1)));
                    nn = nn+1;
                    
                    R(invx,1)=0;
                elseif BC.top(1)==0  % no slip
                    %                  L(invx,invx+3)=Kc/(dy(i)+dy(i+1));
                    %                  L(invx,invx)=-Kc*(1/dy(i)+1/(dy(i)+dy(i+1)));
                    % Multiply whole eq. by 2 so coefficients are closer to
                    % the others
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx+3;
                    Lvv(nn)  = -2*Kc/(dy(i)+dy(i+1));
                    nn = nn+1;
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx;
                    Lvv(nn)  = 2*Kc*(1/dy(i)+1/(dy(i)+dy(i+1)));
                    nn = nn+1;
                    
                    R(invx,1)=2*Kc*(1/dy(i))*BC.top_profile(j);  %G.Ito x-velocity at top boundary;
                end
                
            end
            
            if i==Ny-1 && j<=Nx-1 && j>=2 % X-STOKES: lower boundary not including sides
                if BC.bot(1)==0 % no slip
                    %                  L(invx,invx)=Kc*(1/dy(i)+1/(dy(i-1)+dy(i)));
                    %                  L(invx,invx-3)=-Kc/(dy(i-1)+dy(i));
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx;
                    Lvv(nn)  = 2*Kc*(1/dy(i)+1/(dy(i-1)+dy(i)));
                    nn = nn+1;
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx-3;
                    Lvv(nn)  = -2*Kc/(dy(i-1)+dy(i));
                    nn = nn+1;
                    
                    R(invx,1)=2*Kc*(1/dy(i))*BC.bot_profile(j);
                    
                elseif BC.bot(1)==1 % free slip
                    %                  L(invx,invx)=Kc/(0.5*(dy(i-1)+dy(i)));
                    %                  L(invx,invx-3)=-Kc/(0.5*(dy(i-1)+dy(i)));
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx;
                    Lvv(nn)  = Kc/(0.5*(dy(i-1)+dy(i)));
                    nn = nn+1;
                    
                    Lii(nn) = invx;
                    Ljj(nn) = invx-3;
                    Lvv(nn)  = -Kc/(0.5*(dy(i-1)+dy(i)));
                    nn = nn+1;
                    
                    R(invx,1)=0;
                end
                
                
                
            end
            
            if i==Ny  % lower boundary, ghost vx=0
                %              L(invx,invx)=Kb;
                
                Lii(nn) = invx;
                Ljj(nn) = invx;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                R(invx,1)=0;
            end
            
            
            
        else

            % X-STOKES: internal nodes
            
            % X-STOKES: coeffs for vx
            
            %         L(invx,invx-3*Ny)=4*etan(i+1,j)/(dx(j-1)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invx-3*Ny;
            Lvv(nn)  = 4*etan(i+1,j)/(dx(j-1)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invx,invx)=(-4/(dx(j-1)+dx(j)))*(etan(i+1,j+1)/dx(j)+etan(i+1,j)/dx(j-1))-(2/dy(i))*(etas(i+1,j)/(dy(i)+dy(i+1))+etas(i,j)/(dy(i-1)+dy(i))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invx;
            Lvv(nn)  = (-4/(dx(j-1)+dx(j)))*(etan(i+1,j+1)/dx(j)+etan(i+1,j)/dx(j-1))-(2/dy(i))*(etas(i+1,j)/(dy(i)+dy(i+1))+etas(i,j)/(dy(i-1)+dy(i)));
            nn = nn+1;
            
            %         L(invx,invx+3*Ny)=4*etan(i+1,j+1)/(dx(j)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invx+3*Ny;
            Lvv(nn)  = 4*etan(i+1,j+1)/(dx(j)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invx,invx-3)=2*etas(i,j)/(dy(i)*(dy(i-1)+dy(i))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invx-3;
            Lvv(nn)  = 2*etas(i,j)/(dy(i)*(dy(i-1)+dy(i)));
            nn = nn+1;
            
            %         L(invx,invx+3)=2*etas(i+1,j)/(dy(i)*(dy(i)+dy(i+1))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invx+3;
            Lvv(nn)  = 2*etas(i+1,j)/(dy(i)*(dy(i)+dy(i+1)));
            nn = nn+1;
            
            % X-STOKES: coeffs for vy
            
            %         L(invx,invy+3)=2*etas(i+1,j)/(dy(i)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invy+3;
            Lvv(nn)  = 2*etas(i+1,j)/(dy(i)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invx,invy+3-3*Ny)=-2*etas(i+1,j)/(dy(i)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invy+3-3*Ny;
            Lvv(nn)  = -2*etas(i+1,j)/(dy(i)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invx,invy)=-2*etas(i,j)/(dy(i)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invy;
            Lvv(nn)  = -2*etas(i,j)/(dy(i)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invx,invy-3*Ny)=2*etas(i,j)/(dy(i)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invx;
            Ljj(nn) = invy-3*Ny;
            Lvv(nn)  = 2*etas(i,j)/(dy(i)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            % X-STOKES: coeffs for Pp
            
            %         L(invx,inp+3*Ny+3)=-2*Kc/(dx(j-1)+dx(j)); %
            
            Lii(nn) = invx;
            Ljj(nn) = inp+3*Ny+3;
            Lvv(nn)  = -2*Kc/(dx(j-1)+dx(j));
            nn = nn+1;
            
            %         L(invx,inp+3)=2*Kc/(dx(j-1)+dx(j)); %
            
            Lii(nn) = invx;
            Ljj(nn) = inp+3;
            Lvv(nn)  = 2*Kc/(dx(j-1)+dx(j));
            nn = nn+1;
            
            % RHS
            R(invx,1)=-gx*(rho(i,j)+rho(i+1,j))/2 - ...
                2.*(srhs_xx(i+1,j+1)-srhs_xx(i+1,j))/(dx(j)+dx(j-1)) - ...
                (srhs_xy(i+1,j)-srhs_xy(i,j))/dy(i);
            
        end
        
        
        
        
        %%%%%%% Y-STOKES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        if (j==1 && i<=Ny-1 && i>=2) || (j==Nx-1 && i<=Ny-1 && i>=2) || (j==Nx) || (i==1 && j<=Nx-1) || (i==Ny && j<=Nx-1)
                   
            if j==1 && i>=2 && i<=Ny-1  % left boundary without edges
                if BC.left(1)==1 % free slip
                    %              L(invy,invy)=-2*Kc/(dx(j)+dx(j+1));
                    %              L(invy,invy+3*Ny)=2*Kc/(dx(j)+dx(j+1));
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn)  = -2*Kc/(dx(j)+dx(j+1));
                    nn = nn+1;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy+3*Ny;
                    Lvv(nn)  = 2*Kc/(dx(j)+dx(j+1));
                    nn = nn+1;
                    
                    R(invy,1)=0;
                    
                elseif BC.left(1)==0 % no slip
                    %              L(invy,invy)=-2*Kc*(1/(dx(j)+dx(j+1))+1/dx(j));
                    %              L(invy,invy+3*Ny)=2*Kc/(dx(j+1)+dx(j));
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn)  = -2*Kc*(1/(dx(j)+dx(j+1))+1/dx(j));
                    nn = nn+1;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy+3*Ny;
                    Lvv(nn)  = 2*Kc/(dx(j+1)+dx(j));
                    nn = nn+1;
                    
                    R(invy,1)=0;
                    
                end
            end
            if j==Nx-1 && i<=Ny-1 && i>=2  % right boundary without edges 
                if BC.right(1)==1 % free slip
                    %              L(invy,invy)=Kc/(0.5*(dx(j)+dx(j-1)));
                    %              L(invy,invy-3*Ny)=-Kc/(0.5*(dx(j)+dx(j-1)));
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn)  = Kc/(0.5*(dx(j)+dx(j-1)));
                    nn = nn+1;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy-3*Ny;
                    Lvv(nn)  = -Kc/(0.5*(dx(j)+dx(j-1)));
                    nn = nn+1;
                    
                    R(invy,1)=0;
                    
                elseif BC.right(1)==0 % no slip
                    %              L(invy,invy)=Kc*(1/(dx(j)+dx(j-1))+1/dx(j));
                    %              L(invy,invy-3*Ny)=-Kc/(dx(j)+dx(j-1));
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn)  = Kc*(1/(dx(j)+dx(j-1))+1/dx(j));
                    nn = nn+1;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy-3*Ny;
                    Lvv(nn) = -Kc/(dx(j)+dx(j-1));
                    nn = nn+1;
                    
                    R(invy,1)=0;
                    
                end
                
            end
            if j==Nx  % right boundary, ghost vy=0
                %              L(invy,invy)=Kb;
                
                Lii(nn) = invy;
                Ljj(nn) = invy;
                Lvv(nn)  = Kb;
                nn = nn+1;
                
                R(invy,1)=0;
                
            end
            
            if i==1 && j<=Nx-1  % upper boundary
                if BC.top(2)==0
                    %              L(invy,invy)=Kb;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn)  = Kb;
                    nn = nn+1;
                    
                    R(invy,1)=Kb*BC.top(3);
                elseif BC.top(2)==3
                    %               open top dvy/dy=0, G.Ito
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn) = Kb;
                    nn      = nn+1;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy+3;
                    Lvv(nn) = -Kb;
                    nn      = nn+1;
                    
                    R(invy,1)=0;
                    
                end
            end
            
            if i==Ny && j<=Nx-1  % lower boundary without edges
                
                if BC.bot(2)==0
                    %              L(invy,invy)=Kb;
                    
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn)  = Kb;
                    nn = nn+1;
                    
                    R(invy,1)=Kb*BC.bot(3);
                    
                elseif BC.bot(2) == 2 %%  external boundary
                    
                    dL=300e3;
                                        
                    Lii(nn) = invy;
                    Ljj(nn) = invy;
                    Lvv(nn) = Kc*(1/dL+1/dy(i-1));
                    nn = nn+1;
                    
                    Lii(nn) = invy; 
                    Ljj(nn) = invy-3;
                    Lvv(nn) = Kc*(-1/dy(i-1));
                    nn = nn+1;

                    R(invy, 1) = 0;
                end
                    
            end
            
            
        else
            
            
            
            % internal nodes
            
            % coeffs for vy
            %         L(invy,invy-3*Ny)=2*etas(i,j)/(dx(j)*(dx(j-1)+dx(j))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invy-3*Ny;
            Lvv(nn)  = 2*etas(i,j)/(dx(j)*(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invy,invy)=(-4/(dy(i-1)+dy(i)))*(etan(i+1,j+1)/dy(i)+etan(i,j+1)/dy(i-1))-(2/dx(j))*(etas(i,j+1)/(dx(j)+dx(j+1))+etas(i,j)/(dx(j-1)+dx(j))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invy;
            Lvv(nn)  = (-4/(dy(i-1)+dy(i)))*(etan(i+1,j+1)/dy(i)+etan(i,j+1)/dy(i-1))-(2/dx(j))*(etas(i,j+1)/(dx(j)+dx(j+1))+etas(i,j)/(dx(j-1)+dx(j)));
            nn = nn+1;
            
            %         L(invy,invy+3*Ny)=2*etas(i,j+1)/(dx(j)*(dx(j)+dx(j+1))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invy+3*Ny;
            Lvv(nn)  = 2*etas(i,j+1)/(dx(j)*(dx(j)+dx(j+1)));
            nn = nn+1;
            
            %         L(invy,invy-3)=4*etan(i,j+1)/(dy(i-1)*(dy(i-1)+dy(i))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invy-3;
            Lvv(nn)  = 4*etan(i,j+1)/(dy(i-1)*(dy(i-1)+dy(i)));
            nn  = nn+1;
            
            %         L(invy,invy+3)=4*etan(i+1,j+1)/(dy(i)*(dy(i-1)+dy(i))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invy+3;
            Lvv(nn)  = 4*etan(i+1,j+1)/(dy(i)*(dy(i-1)+dy(i)));
            nn = nn+1;
            
            % coeffs for vx
            
            %         L(invy,invx+3*Ny)=2*etas(i,j+1)/(dx(j)*(dy(i-1)+dy(i)));  %
            
            Lii(nn) = invy;
            Ljj(nn) = invx+3*Ny;
            Lvv(nn)  = 2*etas(i,j+1)/(dx(j)*(dy(i-1)+dy(i)));
            nn = nn+1;
            
            %         L(invy,invx+3*Ny-3)=-2*etas(i,j+1)/(dx(j)*(dy(i-1)+dy(i))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invx+3*Ny-3;
            Lvv(nn)  = -2*etas(i,j+1)/(dx(j)*(dy(i-1)+dy(i)));
            nn = nn+1;
            
            
            %         L(invy,invx)=-2*etas(i,j)/(dx(j)*(dy(i-1)+dy(i))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invx;
            Lvv(nn)  = -2*etas(i,j)/(dx(j)*(dy(i-1)+dy(i)));
            nn = nn+1;
            
            %         L(invy,invx-3)=2*etas(i,j)/(dx(j)*(dy(i-1)+dy(i))); %
            
            Lii(nn) = invy;
            Ljj(nn) = invx-3;
            Lvv(nn)  = 2*etas(i,j)/(dx(j)*(dy(i-1)+dy(i)));
            nn=nn+1;
            
            % coeffs for Pp
            %         L(invy,inp+3*Ny+3)=-2*Kc/(dy(i-1)+dy(i)); %
            
            Lii(nn) = invy;
            Ljj(nn) = inp+3*Ny+3;
            Lvv(nn)  = -2*Kc/(dy(i-1)+dy(i));
            nn=nn+1;
            
            %         L(invy,inp+3*Ny)=2*Kc/(dy(i-1)+dy(i)); %
            
            Lii(nn) = invy;
            Ljj(nn) = inp+3*Ny;
            Lvv(nn)  = 2*Kc/(dy(i-1)+dy(i));
            nn=nn+1;
            
            % RHS (note: sxx = -syy because deviatoric)
            R(invy,1)=-gy*(rho(i,j)+rho(i,j+1))/2 + ...
                2.*(srhs_xx(i+1,j+1)-srhs_xx(i,j+1))/(dy(i)+dy(i-1)) - ...
                (srhs_xy(i,j+1)-srhs_xy(i,j))/dx(j); %
            
        end
        

        
    end % for j
end % for i

nn = nn-1;
%%% Build Sparse Matrix

Li = Lii(1:nn);
Lj = Ljj(1:nn);
Lv = Lvv(1:nn);


L = sparse(Li,Lj,Lv);
% SOLVE LINEAR SYSTEM now down in SiStER_flow_solve G.Ito
% Matlab direct solver
% tic
% S=L\R;
% toc
% 





end

                                                                                                                                                                                                                                                                                                                            SiStER_BC_velocity_profiles.m                                                                       0000640 0002626 0002066 00000003436 13122332436 015225  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   %% ===========================================================================
% SiStER_tangential_velocity_BCs
% Taper tangential velocity BC's so that they meet the velocity conditons
% on the adjoining boundaries gradually.  
%
% This was originally written for wedge models to avoid exporbitant pressure 
% build-up against the left side backstop do to the tangential motion on the base
% G.Ito 7/16
%% ===========================================================================


%-----------------------------------------------------------------------------
% If the BOTTOM is no slip
%-----------------------------------------------------------------------------
if (BC.bot(1)==0)
	if (length(BC.bot)<4); %if a tangential velocity is not defined then its zero
		BC.bot(4)=0;
	end

	BC.bot_profile=BC.bot(4)*ones(Nx,1); 

	if (BC.left(2)==0) %if left side is fixed
		jj=find(x<BC.bot_taperx(1));
		BC.bot_profile(jj)=x(jj).*(BC.bot(4)-BC.left(3))/(BC.bot_taperx(1));
	end

	if (BC.right(2)==0) %if right side is fixed
		jj=find(x>BC.bot_taperx(2));
		BC.bot_profile(jj)=BC.bot(4)+(x(jj)-BC.bot_taperx(2)).*(BC.right(3)-BC.bot(4))/(x(end)-BC.bot_taperx(2));
	end
end

%-----------------------------------------------------------------------------
% If the Bottom is no slip
%-----------------------------------------------------------------------------
if (BC.top(1)==0)
	if (length(BC.top)<4)
		BC.top(4)=0;
	end

	BC.top_profile=BC.top(4)*ones(Nx,1);

    if (BC.left(2)==0)  %if left side is fixed
    	jj=find(x<BC.top_taperx(1))
		BC.top_profile(jj)=BC.left(3) + x(jj).*(BC.top(4)-BC.left(3))/(BC.top_taperx(1));
	end
	if (BC.right(2)==0) %if right side is fixed
		jj=find(x>BC.top_taperx(2))
		BC.top_profile(jj)=BC.top(4) + (x(jj)-BC.top_taperx(2)).*(BC.right(3)-BC.top(4))/(x(end)-BC.top_taperx(2));
	end
end
                                                                                                                                                                                                                                  SiStER_flow_solve.m                                                                                 0000640 0002626 0002066 00000011425 13275723643 013311  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   %==========================================================================
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
        epsII_n=sqrt(0.5*EXX.^2  + 0.5.*(Rn-EXX).^2  + EXY_n.^2 + 0.25.*Rn.^2);  
        epsII_s=sqrt(0.5*EXX_s.^2 +0.5*(Rs-EXX_s).^2 + EXY.^2 + 0.25*Rs.^2);
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
            dum=polyfit([1:nfit]',outit(end-nfit+1:end,1),1);
            conv_rate=abs(dum(1))./mean(outit(end-nfit+1:end,1));
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

      
                                                                                                                                                                                                                                           SiStER_get_cohesion_psim.m                                                                          0000640 0002626 0002066 00000002121 13275245012 014606  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [cohes,psim]=SiStER_marker_cohesion_psim(im,ep,MAT)
% [cohes]=SiStER_get_cohesion(im,ep,MAT)
% compute cohesion on markers based on ep
%G.Ito 8/2016
% Now also computes psim for dilation
% sped up by B. Klein 9/2016
%Added ep_startweakening 4/17

cohes = zeros(size(im));
psim=zeros(size(im));

types = unique(im);
for i = 1:length(types)
    imInd = im == types(i);
    Cmax=MAT(types(i)).Cmax;
    Cmin=MAT(types(i)).Cmin;
    epscrit=MAT(types(i)).ecrit;
    ep1=MAT(types(i)).ep_startweakening; 
    %ep(imInd)=min(ep(imInd),epscrit);
    %disp(['epscrit=',num2str(epscrit)]);
    % get cohesion
    cohes(imInd)=min( max( Cmax+(Cmin-Cmax).*(ep(imInd)-ep1)./(epscrit-ep1),Cmin),Cmax);
    if (isfield(MAT,'psi')~=0)
        ii=find(ep(imInd)>ep1);
        psim(imInd(ii))=MAT(types(i)).psi.*max(1-(ep(imInd(ii))-ep1)./(epscrit-ep1),0);
    else
        psim=0;
    end
end

return

%% OLD VERSION the bracket approach is really slow
% 
% Cmax=[MAT(im).Cmax];
% Cmin=[MAT(im).Cmin];
% epscrit=[MAT(im).ecrit];
% 
% % get cohesion
% cohes=max(Cmax+(Cmin-Cmax).*ep./epscrit,Cmin);
% 
% return


                                                                                                                                                                                                                                                                                                                                                                                                                                               SiStER_get_density.m                                                                                0000640 0002626 0002066 00000001016 13122415241 013422  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [rhom]=SiStER_get_density(im,Tm,MAT)
% [rho]=SiStER_get_density(im,Tm,MAT)
% obtain density from temperature and material identity
% 
% edited by B. Klein, 9/20/2016 to remove struct indexing concatenating

T0=0;
rhom = zeros(size(im));

types = unique(im);
for i = 1:length(types)
    logical = im == types(i);
    rho0 = MAT(types(i)).rho0;
    alpha = MAT(types(i)).alpha;
    rhom(logical) = rho0.*(1-alpha.*(Tm(logical)-T0));
end


% OLD WAY (slow)
%T0=0;
%rhom = [MAT(im).rho0].*(1-[MAT(im).alpha].*(Tm-T0));


                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  SiStER_get_ductile_rheology.m                                                                       0000640 0002626 0002066 00000004041 13122417331 015307  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [eta]=SiStER_get_ductile_rheology(MAT,PARAMS,T,epsII,phase_node)
% 
% computes ductile rheology given phase numbers around a node or cell 
% (phase_node is either phase_n or phase_s)
% Potential shortcoming:  it assumes only <+ TWO phases that are consecutively
% numbered surround each node
% G.Ito 8/2016 % JAO 9/2016 fixed parentheses in exp(Ea/(nRT))
%-------------------------------------------------------------------------
         
[Ny,Nx]=size(epsII);

ph1=floor(phase_node);
ph2=ph1;

if (sum(ceil(phase_node)>ph1+1));
    disp('problems in SiStER_get_ductile_rheology');
    halt
end
pre_diff=cell2mat({MAT(1:end).pre_diff})'; 
ndiff=cell2mat({MAT(1:end).ndiff})';
Ediff=cell2mat({MAT(1:end).Ediff})';
pre_disc=cell2mat({MAT(1:end).pre_disc})'; 
ndisc=cell2mat({MAT(1:end).ndisc})';
Edisc=cell2mat({MAT(1:end).Edisc})';

eta_diff1=pre_diff(ph1).^(-1./ndiff(ph1)).*epsII.^((1-ndiff(ph1))./ndiff(ph1)).*...
          exp(Ediff(ph1)./(ndiff(ph1).*PARAMS.R.*(T+273.15)));
eta_disc1=pre_disc(ph1).^(-1./ndisc(ph1)).*epsII.^((1-ndisc(ph1))./ndisc(ph1)).*...
          exp(Edisc(ph1)./(ndisc(ph1).*PARAMS.R.*(T+273.15)));

eta_diff2=zeros(Ny,Nx);
eta_disc2=zeros(Ny,Nx);
f1=ones(Ny,Nx);

ii=find(phase_node > ph1);
if (length(ii)>1);
ph2(ii)=ph1(ii)+1;
eta_diff2(ii)=pre_diff(ph2(ii)).^(-1./ndiff(ph2(ii))).*epsII(ii).^((1-ndiff(ph2(ii)))./ndiff(ph2(ii))).*...
    exp(Ediff(ph2(ii))./(ndiff(ph2(ii)).*PARAMS.R.*(T(ii)+273.15)));
eta_disc2(ii)=pre_disc(ph2(ii)).^(-1./ndisc(ph2(ii))).*epsII(ii).^((1-ndisc(ph2(ii)))./ndisc(ph2(ii))).*...
    exp(Edisc(ph2(ii))./(ndisc(ph2(ii)).*PARAMS.R.*(T(ii)+273.15)));
f1(ii)=(ph2(ii)-phase_node(ii));  %fraction of phase 1
end;


eta=f1.*(1./eta_diff1+1./eta_disc1).^-1+(1-f1).*(1./eta_diff2+1./eta_disc2).^-1;
    
% if (0);  %testing the shape of cell2mat
%     YY=[[1:10]' [11:20]'];
%     stuff(1).g=1;
%     stuff(2).g=2;
%     ph=YY;
%     ph(:,1)=1; ph(:,2)=2;
%     check1=cell2mat({stuff(ph).g});
%     ZZ=YY(:)'.*check1;
%     ZZ=reshape(ZZ,10,2);
% end
return
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               SiStER_get_elastic_moduli.m                                                                         0000640 0002626 0002066 00000000504 13122415256 014747  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [Gm]=SiStER_get_elastic_moduli(im,MAT)
% [ys]=SiStER_get_elastic_moduli(im,MAT)
% obtain shear modulus from material properties
% B. Klein 9/16

Gm = zeros(size(im));

types = unique(im);
for i = 1:length(types)
    imInd = im == types(i);
    Gm(imInd) = MAT(types(i)).G;
end


% faster than:
% Gm=[MAT(im).G];



                                                                                                                                                                                            SiStER_get_marker_velocities.m                                                                      0000640 0002626 0002066 00000011106 13122410535 015454  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [VX,VY] = SiStER_get_marker_velocities(quad,icn,jcn,x,y,xm,ym,vx,vy,dx,dy)
%
% [VX,VY] = SiStER_get_marker_velocities(quad,icn,jcn,x,y,xm,ym,vx,vy,dx,dy)
% interpolates velocities at nodes to velocities on markers
%
% First cut - J.A. Olive, 2011-2012, modified by B.Z. Klein 2013-2014

Nx = length(x);
Ny = length(y);
[m,~] = size(vx);

%boolean arrays for interpolation
quad1 = quad == 1;
quad2 = quad == 2;
quad3 = quad == 3;
quad4 = quad == 4;

top  = icn == 1;
bot  = icn == Ny-1;

left = jcn == 1;
right= jcn == Nx-1;

%linear index vector of icn & jcn for vx, vy
IND = sub2ind(size(vx), icn, jcn);

%%%%%%%%%%%%%%%%%%%%%%
%% Vx Interpolation %%
%%%%%%%%%%%%%%%%%%%%%%

Vx_xnodes = [x(jcn); x(jcn+1); x(jcn+1); x(jcn)];


%% case 1, icn == 1 (top)
c1 = top;

Vx_ynodes(:,c1) =  [y(icn(c1))+dy(icn(c1))/2;...
                    y(icn(c1))+dy(icn(c1))/2;...
                    y(icn(c1)+1)+dy(icn(c1)+1)/2;...
                    y(icn(c1)+1)+dy(icn(c1)+1)/2];
              
Vx_nodes(:,c1)  =  [vx(IND(c1));...
                    vx(IND(c1)+m);...
                    vx(IND(c1)+1+m);...
                    vx(IND(c1)+1)];

%% case 2, icn == Ny-1 (bot) && in quad 3 or 4
c2 = bot; %& (quad3 | quad4);

Vx_ynodes(:,c2) = [y(icn(c2)-1)+dy(icn(c2)-1)/2;...
                   y(icn(c2)-1)+dy(icn(c2)-1)/2;...
                   y(icn(c2))+dy(icn(c2))/2;...
                   y(icn(c2))+dy(icn(c2))/2];

Vx_nodes(:,c2) = [vx(IND(c2)-1);...
                  vx(IND(c2)-1+m);...
                  vx(IND(c2)+m);...
                  vx(IND(c2))];

%% case 3 icn ~top & ~bot, quad = 1 | 2
c3 = ~top & ~bot & (quad1 | quad2);

Vx_ynodes(:,c3) = [y(icn(c3))-dy(icn(c3)-1)/2;...
                   y(icn(c3))-dy(icn(c3)-1)/2;...
                   y(icn(c3))+dy(icn(c3))/2;...
                   y(icn(c3))+dy(icn(c3))/2];
             
Vx_nodes(:,c3) = [vx(IND(c3)-1);...
                  vx(IND(c3)-1+m);...
                  vx(IND(c3)+m);...
                  vx(IND(c3))];
            
%% case 4 icn ~top & ~bot & quad = 3|4
c4 = ~top & ~bot & (quad3 | quad4);

Vx_ynodes(:,c4) = [y(icn(c4))+dy(icn(c4))/2;...
                   y(icn(c4))+dy(icn(c4))/2;...
                   y(icn(c4)+1)+dy(icn(c4)+1)/2;...
                   y(icn(c4)+1)+dy(icn(c4)+1)/2];
             
Vx_nodes(:,c4) = [vx(IND(c4));...
                  vx(IND(c4)+m);...
                  vx(IND(c4)+1+m);...
                  vx(IND(c4)+1)];

            
            
VX = SiStER_interp_grid_to_marker_vector(Vx_xnodes, Vx_ynodes, Vx_nodes, xm, ym);



%%%%%%%%%%%%%%%%%%%%%%%
%% VY INTERPOLATION  %%
%%%%%%%%%%%%%%%%%%%%%%%

Vy_ynodes = [y(icn); y(icn); y(icn+1); y(icn+1)];

%% case 1 jcn == 1 (left)
c1 = left;

Vy_xnodes(:,c1) = [x(jcn(c1))+dx(jcn(c1))/2;...
                   x(jcn(c1)+1)+dx(jcn(c1)+1)/2;...
                   x(jcn(c1)+1)+dx(jcn(c1)+1)/2;...
                   x(jcn(c1))+dx(jcn(c1))/2];
           
Vy_nodes(:,c1) = [vy(IND(c1));...
                  vy(IND(c1)+m);...
                  vy(IND(c1)+1+m);...
                  vy(IND(c1)+1)];

%% case 2 jcn == Nx-1 (right) & Quad = 2 | 3
c2 = right;% & (quad2 | quad3);???????????????????????????

Vy_xnodes(:,c2) = [x(jcn(c2)-1)+dx(jcn(c2)-1)/2;...
                   x(jcn(c2))+dx(jcn(c2))/2;...
                   x(jcn(c2))+dx(jcn(c2))/2;...
                   x(jcn(c2)-1)+dx(jcn(c2)-1)/2];
           
Vy_nodes(:,c2) = [vy(IND(c2)-m);...
                  vy(IND(c2));...
                  vy(IND(c2)+1);...
                  vy(IND(c2)+1-m)];
            
 %% case 3 ~left & ~right & quad = 1 | 4
 c3 = ~left & ~right & (quad1 | quad4);
 
 Vy_xnodes(:,c3) = [x(jcn(c3))-dx(jcn(c3)-1)/2;...
                    x(jcn(c3))+dx(jcn(c3))/2;...
                    x(jcn(c3))+dx(jcn(c3))/2;...
                    x(jcn(c3))-dx(jcn(c3)-1)/2];
              
 Vy_nodes(:,c3) = [vy(IND(c3)-m);...
                   vy(IND(c3));...
                   vy(IND(c3)+1);...
                   vy(IND(c3)+1-m)];
            
%% case 4 ~left & ~right & quad = 2 | 3
c4 = ~left & ~right & (quad2 | quad3);

Vy_xnodes(:,c4) = [x(jcn(c4))+dx(jcn(c4))/2;...
                   x(jcn(c4)+1)+dx(jcn(c4)+1)/2;...
                   x(jcn(c4)+1)+dx(jcn(c4)+1)/2;...
                   x(jcn(c4))+dx(jcn(c4))/2];
             
Vy_nodes(:,c4) = [vy(IND(c4));...
                  vy(IND(c4)+m);...
                  vy(IND(c4)+m+1);...
                  vy(IND(c4)+1)];
 

            
VY = SiStER_interp_grid_to_marker_vector(Vy_xnodes,Vy_ynodes,Vy_nodes, xm, ym);            
                                                                                                                                                                                                                                                                                                                                                                                                                                                          SiStER_get_mean_material_prop.m                                                                     0000640 0002626 0002066 00000001707 13122410535 015611  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   %==========================================================================
function [var] = SiStER_get_mean_material_prop(matp,phase_node)
%
% Computes mean value of MAT.{property} at nodes
% e.g., [Mu_s]=SiStER_material_prop(cell2mat({MAT([1:Nphase]).mu}),phase_s);
%
% Shortcoming:  phase_node (=phase_s or phase_n) is phase number determined
% by interpolating im from marker to (shear or normal) node.  Thus this will
% work only if there are two phases adjacent to the nodes and those two 
% phases are consecutively numbered. 
% G.Ito 8/16
%==========================================================================


phase1=floor(phase_node);
var=matp(phase1);

ii=find(phase_node > phase1);
if (length(ii)>1);
    if (sum(ceil(phase_node)>phase1+1));
        disp('problems in SiStER_material_prop_s');
        halt
    end;
   
    f1=((phase1(ii)+1)-phase_node(ii))';  %fraction of phase 1
    var(ii)=f1.*matp(phase1(ii))+(1-f1).*matp(phase1(ii)+1);
end

                                                         SiStER_get_mu.m                                                                                     0000640 0002626 0002066 00000001120 13136445620 012371  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [mum]=SiStER_get_mu(im,ep,MAT)
% [cohes]=SiStER_get_cohesion(im,ep,MAT)
% compute coefficient of friction on markers based on ep
%G.Ito 7/2017
% using B. Klein's efficient use of structures w/markers 9/2016

mum = zeros(size(im));

types = unique(im);
for i = 1:length(types)
    imInd = im == types(i);
    mumax=MAT(types(i)).mu;
    mumin=MAT(types(i)).mumin;
    epscrit=MAT(types(i)).ecrit;
    ep1=MAT(types(i)).ep_startweakening; 

    %weakened friction on each marker
    mum(imInd)=min( max( mumax+(mumin-mumax).*(ep(imInd)-ep1)./(epscrit-ep1),mumin),mumax);
end

return


                                                                                                                                                                                                                                                                                                                                                                                                                                                SiStER_get_rotation_rate.m                                                                          0000640 0002626 0002066 00000006546 13122410535 014633  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [om]=SiStER_get_rotation_rate(vx,vy,dx,dy,BC)
% [EXX EXY]=HIPSTER_get_rotation_rate(vx,vy,dx,dy,BC)
% computes the rotation rate at shear nodes

[Ny, Nx]=size(vx);
om=zeros(size(vy));


for i=1:Ny
    for j=1:Nx   
        
        if i>=2 && i<=Ny-1 && j>=2 && j<=Nx-1
            om(i,j)=0.5*(-2*(vx(i,j)-vx(i-1,j))/(dy(i-1)+dy(i)) + 2*(vy(i,j)-vy(i,j-1))/(dx(j-1)+dx(j)));
        end
        
        if i==1 && j>=2 && j<=Nx-1 % top
            if BC.top(1)==1 % rollers
                dvxdy=0;
            elseif BC.top(1)==0 % no rollers
                dvxdy=2*vx(i,j)/dy(i);
            end         
            om(i,j)= 0.5*(-dvxdy + 2*(vy(i,j)-vy(i,j-1))/(dx(j-1)+dx(j)));
            
        elseif i==Ny && j>=2 && j<=Nx-1 % bottom
            if BC.bot(1)==1 % rollers
                dvxdy=0;
            elseif BC.bot(1)==0 % no rollers
                dvxdy=-2*vx(i-1,j)/dy(i-1);
            end         
            om(i,j)= 0.5*(-dvxdy + 2*(vy(i,j)-vy(i,j-1))/(dx(j-1)+dx(j)));
            
        elseif j==1 && i>=2 && i<=Ny-1
            if BC.left(1)==1 % rollers
                dvydx=0;
            elseif BC.left(1)==0 % no rollers
                dvydx=2*vy(i,j)/dx(j);
            end     
            om(i,j)=0.5*(-2*(vx(i,j)-vx(i-1,j))/(dy(i-1)+dy(i)) + dvydx);
            
        elseif j==Nx && i>=2 && i<=Ny-1
            if BC.right(1)==1 % rollers
                dvydx=0;
            elseif BC.right(1)==0 % no rollers
                dvydx=-2*vy(i,j-1)/dx(j-1);
            end     
            om(i,j)=0.5*(-2*(vx(i,j)-vx(i-1,j))/(dy(i-1)+dy(i)) + dvydx);
            
        end
        
        % corners
        if i==1 && j==1
            
            if BC.top(1)==1 % rollers
                dvxdy=0;
            elseif BC.top(1)==0 % no rollers
                dvxdy=2*vx(i,j)/dy(i);
            end   
            if BC.left(1)==1 % rollers
                dvydx=0;
            elseif BC.left(1)==0 % no rollers
                dvydx=2*vy(i,j)/dx(j);
            end   
            om(i,j)=(-dvxdy + dvydx)/2;
            
        elseif i==1 && j==Ny
            
            if BC.top(1)==1 % rollers
                dvxdy=0;
            elseif BC.top(1)==0 % no rollers
                dvxdy=2*vx(i,j)/dy(i);
            end  
            if BC.right(1)==1 % rollers
                dvydx=0;
            elseif BC.right(1)==0 % no rollers
                dvydx=-2*vy(i,j-1)/dx(j-1);
            end  
            om(i,j)=(-dvxdy + dvydx)/2;
            
        elseif i==Nx && j==Ny
            if BC.right(1)==1 % rollers
                dvydx=0;
            elseif BC.right(1)==0 % no rollers
                dvydx=-2*vy(i,j-1)/dx(j-1);
            end  
            if BC.bot(1)==1 % rollers
                dvxdy=0;
            elseif BC.bot(1)==0 % no rollers
                dvxdy=-2*vx(i-1,j)/dy(i-1);
            end  
            om(i,j)=(-dvxdy + dvydx)/2;
            
        elseif i==Nx && j==1
            if BC.bot(1)==1 % rollers
                dvxdy=0;
            elseif BC.bot(1)==0 % no rollers
                dvxdy=-2*vx(i-1,j)/dy(i-1);
            end  
            if BC.left(1)==1 % rollers
                dvydx=0;
            elseif BC.left(1)==0 % no rollers
                dvydx=2*vy(i,j)/dx(j);
            end  
            om(i,j)=(-dvxdy + dvydx)/2;
            
        end
            
            
        
    end
end                                                                                                                                                                  SiStER_get_strain_rate.m                                                                            0000640 0002626 0002066 00000012304 13275755274 014306  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [EXX,EXY]=SiStER_get_strain_rate(vx,vy,dx,dy,BC)
% [EXX EXY]=SiStER_get_strain_rate(vx,vy,dx,dy,BC)
% computes the deviatoric strain rate: EXX, at normal nodes
% and EXY, at shear nodes, from vx and vy on the eulerian grid
%
% G.Ito 11/12/15 updated calculations at corner nodes to be
%    consistent with the top and bottom boundary conditions, as those
%    (not the side) conditions control the solutions for vx,vy at the
%    points closest to the corners in SiStER_Stokes_solver.m

[Ny, Nx]=size(vx);

EXX=zeros(size(vx));
EXY=zeros(size(vy));

for i=2:Ny
    for j=2:Nx   
        EXX(i,j)=(vx(i-1,j)-vx(i-1,j-1))/dx(j-1);        
    end
end
        

for i=1:Ny
    for j=1:Nx   
        
        if i>=2 && i<=Ny-1 && j>=2 && j<=Nx-1
            EXY(i,j)=0.5*(2*(vx(i,j)-vx(i-1,j))/(dy(i-1)+dy(i)) + 2*(vy(i,j)-vy(i,j-1))/(dx(j-1)+dx(j)));
        end
        
        if i==1 && j>=2 && j<=Nx-1 % top
            if BC.top(1)==1 % rollers
                dvxdy=0;
            elseif BC.top(1)==0 % no rollers
                dvxdy=2*(vx(i,j)-BC.top_profile(j))/dy(i);  %G.Ito tangential velocity on top 5/12/16
            end         
            EXY(i,j)= 0.5*(dvxdy + 2*(vy(i,j)-vy(i,j-1))/(dx(j-1)+dx(j)));
            
        elseif i==Ny && j>=2 && j<=Nx-1 % bottom
            if BC.bot(1)==1 % rollers
                dvxdy=0;
            elseif BC.bot(1)==0 % no rollers
                %dvxdy=-2*vx(i-1,j)/dy(i-1);
                dvxdy=2.*(BC.bot_profile(j)-vx(i-1,j))/dy(i-1);  %G.Ito needed for tangential velocity on base
            end         
            EXY(i,j)= 0.5*(dvxdy + 2*(vy(i,j)-vy(i,j-1))/(dx(j-1)+dx(j)));
            
        elseif j==1 && i>=2 && i<=Ny-1 
            if BC.left(1)==1 % rollers
                dvydx=0;
            elseif BC.left(1)==0 % no rollers
                dvydx=2*vy(i,j)/dx(j);
            end     
            EXY(i,j)=0.5*(2*(vx(i,j)-vx(i-1,j))/(dy(i-1)+dy(i)) + dvydx);
            
        elseif j==Nx && i>=2 && i<=Ny-1
            if BC.right(1)==1 % rollers
                dvydx=0;
            elseif BC.right(1)==0 % no rollers
                dvydx=-2*vy(i,j-1)/dx(j-1);
            end     
            EXY(i,j)=0.5*(2*(vx(i,j)-vx(i-1,j))/(dy(i-1)+dy(i)) + dvydx);
            
        end
        
        % corners
        if i==1 && j==1
            
            if BC.top(1)==1 % rollers
                dvxdy=0;
            elseif BC.top(1)==0 % no rollers
                %dvxdy=2*vx(i,j)/dy(i);
                dvxdy=2*(vx(i,j)-BC.top_profile(j))/dy(i);  %G.Ito tangential velocity on top 5/12/16
            end   
            
            if BC.top(2)==0 % vy imposed on top (Dirichlet)
                dvydx=0;
            elseif BC.top(2)==3 % open top so set to be the same as adjacent node
                dvydx=2*(vy(i,j+1)-vy(i,j))/(dx(j)+dx(j+1));
            else
                disp(['Stopped in get_strain_rate:  EXY for corner elements not coded for BC.top(2)='....
                num2str(BC.top(2))]);
                halt
            end   
            EXY(i,j)=(dvxdy + dvydx)/2;
            
        elseif i==1 && j==Nx 
            
            if BC.top(1)==1 % rollers
                dvxdy=0;
            elseif BC.top(1)==0 % no slip
                %dvxdy=2*vx(i,j)/dy(i);
                dvxdy=2*(vx(i,j)-BC.top_profile(j))/dy(i);  %G.Ito tangential velocity on top
            end  
            if BC.top(2)==0 % vy imposed on top (Dirichlet)
                dvydx=0;
            elseif BC.top(2)==3 % open top so set to be the same as adjacent node
                 dvydx=2*(vy(i,j-1)-vy(i,j-2))./(dx(j-1)+dx(j-2));
            end  
            EXY(i,j)=(dvxdy + dvydx)/2;
            
        elseif i==Ny && j==Nx
            if BC.bot(1)==1 % rollers
                dvxdy=0;
            elseif BC.bot(1)==0 % no rollers
                dvxdy=2.*(BC.bot_profile(j)-vx(i-1,j))/dy(i-1);  %G.Ito 5/11/16 corner bottom same as left side    
            end
            if BC.bot(2)==0 % vy imposed on bottom (Dirichlet)
                dvydx=0;
            elseif BC.bot(2)==3 % open bottom so set to be the same as adjacent node
                dvydx=2*(vy(i,j-1)-vy(i,j-2))./(dx(j-1)+dx(j-2));
            else
                disp(['Stopped in get_strain_rate:  EXY for corner elements not coded for BC.bot(2)='....
                num2str(BC.bot(2))]);
                halt
            end
            EXY(i,j)=(dvxdy + dvydx)/2;
            
        elseif i==Ny && j==1
            if BC.bot(1)==1 % rollers
                dvxdy=0;
            elseif BC.bot(1)==0 % no rollers
                %dvxdy=-2*vx(i-1,j)/dy(i-1);
                dvxdy=2.*(BC.bot_profile(j)-vx(i-1,j))/dy(i-1);  %G.Ito 5/11/16 corner bottom same as left side            
            end  
            if BC.bot(2)==0 % vy imposed on bottom (Dirichlet)% rollers
                dvydx=0;
            elseif BC.bot(2)==3 % open bottom so set to be the same as adjacent node
                dvydx=2*(vy(i,j+1)-vy(i,j))./(dx(j+1)+dx(j));
            end  
            EXY(i,j)=(dvxdy + dvydx)/2;
            
        end
            
            
        
    end
end                                                                                                                                                                                                                                                                                                                                    SiStER_get_thermal_properties.m                                                                     0000640 0002626 0002066 00000000273 13122410535 015660  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [km, cpm]=SiStER_get_thermal_properties(im,MAT)
% [ys]=SiStER_get_thermal_properties(im,MAT)
% obtain thermal conductivity and heat capacity

km=[MAT(im).k];
cpm=[MAT(im).cp];


                                                                                                                                                                                                                                                                                                                                     SiStER_initialize_grid_GI.m                                                                         0000640 0002626 0002066 00000003224 13122410535 014635  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [X,Y,x,y,xc,yc,dx,dy,Nx,Ny] = SiStER_initialize_grid_GI(xsize,ysize,GRID)

%% Coded for both regular and variable grid spacing G.Ito 7/15
%
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

%% The same logic is used for y.


x=0;
y=0;
for i=1:length(GRID.x);
  n=round((GRID.x(i)-x(end))/GRID.dx(i));  %set number of element in zone i
  dd=(GRID.x(i)-x(end))/n;                 %set real grid size
  x=[x, x(end)+(1:n).*dd];
end
n=round((xsize-x(end))./GRID.dx(end));
dd=(xsize-x(end))./n;
x=[x, x(end)+(1:n).*dd];

for i=1:length(GRID.y);
  n=round((GRID.y(i)-y(end))./GRID.dy(i));
  dd=(GRID.y(i)-y(end))./n;
  y=[y, y(end)+(1:n).*dd];
end
n=round((ysize-y(end))./GRID.dy(end));
dd=(ysize-y(end))./n;
y=[y, y(end)+(1:n).*dd];

Nx=length(x);
Ny=length(y);  % they may have changed !

% 2D Grid
[X,Y] = meshgrid(x,y);

% create vectors of grid spacing
dx=diff(x);
dy=diff(y);

% Create vectors for cell centers positions (staggered nodes)
xc=0.5*(x(2:Nx)+x(1:Nx-1));
yc=0.5*(y(2:Ny)+y(1:Ny-1));

disp(['** Nx = ' num2str(Nx) ' Ny = ' num2str(Ny) ' ***']);


                                                                                                                                                                                                                                                                                                                                                                            SiStER_Initialize.m                                                                                 0000640 0002626 0002066 00000014122 13275724433 013226  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   %===========================================================================
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
if (exist('strata_y','var'));                                          %G.Ito
    [stratam]=SiStER_initialize_marker_stratigraphy(strata_y,xm,ym); %G.Ito
else
    stratam=zeros(size(sxxm));                                      %G.Ito
end;                                                                %G.Ito

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

Rn=zeros(Ny,Nx); %This must be initialized for assemble_L_R


%%
%Set minimum number of Picard iterations
if isfield(PARAMS,'Npicard_min')==0;
         PARAMS.Npicard_min=1;
end

% FAULTING SEED
% randomly seeds ep to promote faulting if switch is on (internal)
if (exist('PARAMS.plastic_seed','var') && sum(PARAMS.plastic_seed ~= 'none  ')<6);
SiStER_plastic_seed;
end;

%% -------------------------------------------------------------------------
% SURFACE MARKERS  G.Ito
%-------------------------------------------------------------------------
if exist('isurface','var')
    if isurface==1
        Ntopo=round(xsize/PARAMS.dxsurf)+1;
        topo_x=linspace(0,xsize,Ntopo);
        topo_y=GEOM(1).bot*ones(size(topo_x));
        if (isfield(BC,'bot_xbackstop') && BC.bot_xbackstop(2)>0);
            kk=4;
            dum1=GEOM(kk).ytopleft+(GEOM(kk).yright-GEOM(kk).ytopleft)/(GEOM(kk).xright-GEOM(kk).xtopleft)*(topo_x-GEOM(kk).xtopleft);
            dum1=min([dum1' GEOM(kk).yright.*ones(size(topo_y')) topo_y']');
            topo_y=dum1;
        end;
    end
end

%% -------------------------------------------------------------------------
% initialize dt_m small to keep things elastic & no plasticity at t=1, G.Ito
%-------------------------------------------------------------------------
if (exist('dt_m','var')==0);
    dt_m=1e2;
end
sxxOLD=EXX; sxyOLD=EXY;  %inialize these as zeros G.Ito


%% --------------------------------------------------------------------------
% Set DEFAULT params G.Ito
%--------------------------------------------------------------------------
if (isfield(PARAMS,'plast_option')==0)
    PARAMS.plast_option=2;           %1=harmonic sum of eta_plas and eta_new, 2=limit eta_new to eta_plas in VEP_rheology
end;
if (isfield(PARAMS,'welast')==0)
    PARAMS.welast=0;           %1=harmonic sum of eta_plas and eta_new, 2=limit eta_new to eta_plas in VEP_rheology
end;

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



                                                                                                                                                                                                                                                                                                                                                                                                                                              SiStER_initialize_marker_phases_backstop.m                                                          0000640 0002626 0002066 00000002031 13122410535 020036  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [im] = SiStER_initialize_marker_phases(Nphase,GEOM,xm,ym,PARAMS,xsize)

% assign material identity on markers
im=zeros(size(xm));


for kk = 1:3
    
    if GEOM(kk).type==1 % layer
        
        im(ym>=GEOM(kk).top & ym<GEOM(kk).bot)=kk;
        
    
    elseif GEOM(kk).type==2 % circular inclusion
        
        rm=sqrt((xm-GEOM(kk).x0).^2 + (ym-GEOM(kk).y0).^2);
        im(rm<=GEOM(kk).rad)=kk;
    end;
    
    
end

if (Nphase==4);
    kk=Nphase;
    ytop=GEOM(kk).ytopleft+(GEOM(kk).yright-GEOM(kk).ytopleft)/(GEOM(kk).xright-GEOM(kk).xtopleft)*(xm-GEOM(kk).xtopleft);
    ytop=min([ytop' GEOM(kk).yright.*ones(size(ym'))]');
    im(ym>=ytop)=kk;    
end

if PARAMS.ridge == 1
    im((ym <= (GEOM(2).bot+((GEOM(2).bot-GEOM(2).top)/PARAMS.Ldouble)*(xm-xsize/2))) ...
       &ym >= GEOM(kk).top ...
       &xm >  xsize/2) = 2;
    im((ym <= (GEOM(2).bot+((GEOM(2).bot-GEOM(2).top)/PARAMS.Ldouble)*(xsize/2-xm))) ...
       &ym >= GEOM(kk).top ...
       &xm <  xsize/2) = 2;
   
end

end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       SiStER_initialize_marker_phases.m                                                                   0000640 0002626 0002066 00000001705 13122410535 016157  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [im] = SiStER_initialize_marker_phases(Nphase,GEOM,xm,ym,PARAMS,xsize)

% assign material identity on markers
im=zeros(size(xm));


for kk = 1:Nphase
    
    if GEOM(kk).type==1 % layer
        
        im(ym>=GEOM(kk).top & ym<GEOM(kk).bot)=kk;
        
    
    elseif GEOM(kk).type==2 % circular inclusion
        rm=sqrt((xm-GEOM(kk).x0).^2 + (ym-GEOM(kk).y0).^2);
        im(rm<=GEOM(kk).rad)=kk;
        
    elseif GEOM(kk).type==3 %rectangular inclusion
        im(GEOM(kk).left<=xm & xm<=GEOM(kk).right & GEOM(kk).top <= ym & ym<=GEOM(kk).bot)=kk;
        
    end
    
end
        
if PARAMS.ridge == 1
    im((ym <= (GEOM(2).bot+((GEOM(2).bot-GEOM(2).top)/PARAMS.Ldouble)*(xm-xsize/2))) ...
       &ym >= GEOM(kk).top ...
       &xm >  xsize/2) = 2;
    im((ym <= (GEOM(2).bot+((GEOM(2).bot-GEOM(2).top)/PARAMS.Ldouble)*(xsize/2-xm))) ...
       &ym >= GEOM(kk).top ...
       &xm <  xsize/2) = 2;
   
end

end                                                           SiStER_initialize_marker_positions.m                                                                0000640 0002626 0002066 00000001264 13122410535 016723  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [xm, ym] = SiStER_initialize_marker_positions(xsize,ysize,dx,dy,Mquad)
% [xm ym dysub] = HIPSTER_initialize_markers(xsize,ysize,dx,dy,mark_per_quadrant)
% assigns markers their inital position



% smallest quadrant sizes
mdx=min(dx)/2;
mdy=min(dy)/2;

% creating a regular grid with step = that of the smallest quadrant
xx=0:mdx:xsize;
yy=0:mdy:ysize;

nxx=length(xx);
nyy=length(yy);


midx=1;
for i=1:nyy-1
    for j=1:nxx-1
            
        [xmtemp, ymtemp]=SiStER_seed_markers_uniformly(xx(j),yy(i),mdx,mdy,Mquad);
        xm(midx:midx+Mquad-1)=xmtemp;
        ym(midx:midx+Mquad-1)=ymtemp;
        midx=midx+Mquad;    
        
    end
end


                                                                                                                                                                                                                                                                                                                                            SiStER_initialize_marker_stratigraphy.m                                                             0000640 0002626 0002066 00000000707 13122410535 017416  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [stratam] = SiStER_initialize_marker_strata(strata_y,xm,ym)

% assign passive stratigraphic horizons
% these do not change material properties, just used to show deformation
% G.Ito 7/15
%-------------------------------------------------------------------------
stratam=zeros(size(xm));


for kk = 1:length(strata_y)-1
    stratam(ym>=strata_y(kk) & ym< strata_y(kk+1) )=kk;
end

stratam(ym>=strata_y(end)) = length(strata_y);

end                                                         SiStER_interp_grid_to_marker_vector.m                                                               0000640 0002626 0002066 00000001077 13122410535 017047  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [param] = SiStER_interp_grid_to_marker_vector(xnodes,ynodes,paramnodes,xm,ym)
% [param] = HIPSTER_interp_grid_to_marker_vector(xnodes,ynodes,paramnodes,xm,ym)
% interpolates a grid value to the entire vector of markers,
% used by SiStER_get_marker_velocities
% B.Z. Klein, 2013

dxm=(xm-xnodes(1,:));
dym=(ym-ynodes(1,:));
dx=abs(xnodes(2,:)-xnodes(1,:));
dy=abs(ynodes(3,:)-ynodes(2,:));

w1=(1-dxm./dx).*(1-dym./dy);
w2=(dxm./dx).*(1-dym./dy);
w4=(dym./dy).*(1-dxm./dx);
w3=(dxm./dx).*(dym./dy);
wnodes=[w1; w2; w3; w4];
    
param=sum(paramnodes.*wnodes, 1);

end                                                                                                                                                                                                                                                                                                                                                                                                                                                                 SiStER_interp_markers_to_normal_nodes.m                                                             0000640 0002626 0002066 00000002005 13122410535 017373  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,varargin)

% [n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,varargin)
% interpolates parameters (in the order of input) from markers to normal nodes
%
% First cut - J.A. Olive, March 2011
% Modified by E. Mittelstaedt, April 2011 to allow multiple inputs  
% Modified by B.Z. Klein, Spring 2014 for speedup

Nx=length(x);
Ny=length(y);
dx=diff(x);
dy=diff(y);

% check for number of properties to interpolate
numV = size(varargin,2);
n2interp(1:numV) = struct('data',zeros(Ny,Nx));

INDEX = sub2ind([Ny-1, Nx-1], icn, jcn);

AcCell = bsxfun(@times, dy', dx);

xN = x(1:Nx-1) + dx/2;
yN = y(1:Ny-1) + dy/2;
[XN, YN] = meshgrid(xN, yN);

AMvec = abs((xm - XN(INDEX)).*(ym - YN(INDEX)));
WMvec = (AcCell(INDEX) - AMvec)./AcCell(INDEX);

w = accumarray([icn' jcn'], WMvec', [], @sum);

for vn = 1:numV
    VecData = varargin{vn}.*WMvec;
    n2interp(vn).data(2:Ny,2:Nx) = accumarray([icn' jcn'], VecData', [], @sum)./w;
end

end




                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           SiStER_interp_markers_to_shear_nodes.m                                                              0000640 0002626 0002066 00000017611 13122410535 017216  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [n2interp] = SiStER_interp_markers_to_shear_nodes(xm,ym,icn,jcn,quad,x,y,varargin)

% First cut - J.A. Olive, March 2011
% Modified by E. Mittelstaedt, April 2011, to allow variable inputs.  
% Modified by B.Z. Klein, Spring 2014, for speedup
% Modified by B.Z. Klein, Summer 2014, for further speedup (vectorized)


Nx=length(x);
Ny=length(y);
dx=diff(x);
dy=diff(y);

% MITTELSTAEDT - check for number of properties to interpolate
numV = size(varargin,2);


% MITTELSTAEDT
% establish interpolants matrices
n2interp = repmat(struct('data', zeros(Ny,Nx)), 1, numV);

weight_cells=0; % weight cells surrounding one node based on their relative areas ? 1=yes - 0=no (Gerya's way)

JCN = interp1(x, 1:length(x), xm, 'nearest', 'extrap'); %% these are like the jcn and icn elsewhere, except the nodes are centered instead of upper left.
ICN = interp1(y, 1:length(y), ym, 'nearest', 'extrap'); %% this makes a lot of the indexing much simpler below.


%% Interior Cells

center = jcn>1 & jcn<Nx & icn>1 & icn<Ny;
shiftLeft = jcn<Nx-1 & icn>1 & icn<Ny;
shiftUp = jcn>1 & jcn<Nx & icn<Ny-1;
shiftBoth = jcn<Nx-1 & icn<Ny-1;


cell1 = center & ((xm-x(JCN)) > 0) & ((ym - y(ICN)) > 0);  %% these are logical arrays that index the original quadrants
cell2 = shiftLeft & ((xm-x(JCN)) < 0) & ((ym - y(ICN)) > 0);
cell3 = shiftBoth & ((xm-x(JCN)) < 0) & ((ym - y(ICN)) < 0);
cell4 = shiftUp & ((xm-x(JCN)) > 0) & ((ym - y(ICN)) < 0);

%%% WEIGHTING (equal for now because that is what I'm running)

wc1 = 0.25;
wc2 = 0.25;
wc3 = 0.25;
wc4 = 0.25;


% cell 1 (i,j,1)

dxm = xm(cell1) - x(JCN(cell1));
dym = ym(cell1) - y(ICN(cell1));
ddx = dx(JCN(cell1));
ddy = dy(ICN(cell1));

wm1 = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w1 = accumarray([ICN(cell1)', JCN(cell1)'], wm1);

% cell 2 (i, j-1, 2)

dxm = xm(cell2) - x(JCN(cell2));
dym = ym(cell2) - y(ICN(cell2));
ddx = dx(JCN(cell2)-1);
ddy = dy(ICN(cell2));

wm2 = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w2 = accumarray([ICN(cell2)', JCN(cell2)'], wm2);

% cell 3 (i-1, j-1, 3)

dxm = xm(cell3) - x(JCN(cell3));
dym = ym(cell3) - y(ICN(cell3));
ddx = dx(JCN(cell3)-1);
ddy = dy(ICN(cell3)-1);

wm3 = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w3 = accumarray([ICN(cell3)', JCN(cell3)'], wm3);

% cell 4 (i-1, j, 4)

dxm = xm(cell4) - x(JCN(cell4));
dym = ym(cell4) - y(ICN(cell4));
ddx = dx(JCN(cell4));
ddy = dy(ICN(cell4)-1);

wm4 = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w4 = accumarray([ICN(cell4)', JCN(cell4)'], wm4);

%loop over material properties to interpolate

for vn = 1:numV
    n2interp(vn).data = (wc1*accumarray([ICN(cell1)', JCN(cell1)'], varargin{vn}(cell1).*wm1)./w1 + ...
        wc2*accumarray([ICN(cell2)', JCN(cell2)'], varargin{vn}(cell2).*wm2)./w2 + ...
        wc3*accumarray([ICN(cell3)', JCN(cell3)'], varargin{vn}(cell3).*wm3)./w3 + ...
        wc4*accumarray([ICN(cell4)', JCN(cell4)'], varargin{vn}(cell4).*wm4)./w4)./...
        (wc1+wc2+wc4+wc4);
end



%% EDGES

%%% top edge

topEdge = jcn>1 & jcn<Nx & icn==1;
shifted = jcn<Nx-1 & icn==1;

% cell 1

cell1 = shifted & quad==2;

ddx = dx(JCN(cell1)-1);
ddy = dy(1);
dxm = xm(cell1) - x(JCN(cell1));
dym = ym(cell1) - y(ICN(cell1));
wm1 = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w1 = accumarray([ICN(cell1)', JCN(cell1)'], wm1);

% cell 2 

cell2 = topEdge & quad==1;

ddx = dx(JCN(cell2));
ddy = dy(1);
dxm = xm(cell2) - x(JCN(cell2));
dym = ym(cell2) - y(ICN(cell2));
wm2 = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w2  = accumarray([ICN(cell2)', JCN(cell2)'], wm2);

%loop over material properties to interpolate

for vn = 1:numV
    temp = (wc1*accumarray([ICN(cell1)', JCN(cell1)'], varargin{vn}(cell1).*wm1)./w1 + ...
        wc2*accumarray([ICN(cell2)', JCN(cell2)'], varargin{vn}(cell2).*wm2)./w2)/...
        (wc1+wc2);
    n2interp(vn).data(1,2:end) = temp(2:end);
end

clear w1 w2

%%% bottom edge

bottomEdge = jcn>1 & jcn<Nx & icn==Ny-1;
shifted =    jcn<Nx-1       & icn==Ny-1;

% cell 1

cell1 = shifted & quad==3;

ddx = dx(JCN(cell1)-1);
ddy = dy(Ny-1);
dxm = xm(cell1) - x(JCN(cell1));
dym = ym(cell1) - y(end-1);
wm1 = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w1 = accumarray([ones(sum(cell1),1), JCN(cell1)'], wm1);

% cell 2

cell2 = bottomEdge & quad==4;

ddx = dx(JCN(cell2));
ddy = dy(Ny-1);
dxm = xm(cell2) - x(JCN(cell2));
dym = ym(cell2) - y(end-1);
wm2 = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w2  = accumarray([ones(sum(cell2),1), JCN(cell2)'], wm2);

%loop over material properties to interpolate

for vn = 1:numV
    temp = (wc1*accumarray([ones(sum(cell1),1), JCN(cell1)'], varargin{vn}(cell1).*wm1)./w1 + ...
        wc2*accumarray([ones(sum(cell2),1), JCN(cell2)'], varargin{vn}(cell2).*wm2)./w2)/...
        (wc1+wc2);
    n2interp(vn).data(Ny,2:end) = temp(2:end);
end

%%% left edge

leftEdge = jcn==1 & icn>1 & icn<Ny;
shifted  = jcn==1 & icn<Ny-1;

% cell 1

cell1 = shifted & quad==4;

ddx = dx(1);
ddy = dy(ICN(cell1)-1);
dxm = xm(cell1) - x(1);
dym = ym(cell1) - y(ICN(cell1));
wm1 = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w1 = accumarray([ICN(cell1)', ones(sum(cell1),1)], wm1);

% cell 2

cell2 = leftEdge & quad==1;

ddx = dx(1);
ddy = dy(ICN(cell2));
dxm = xm(cell2) - x(1);
dym = ym(cell2) - y(ICN(cell2));
wm2 = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w2 = accumarray([ICN(cell2)', ones(sum(cell2),1)], wm2);

%loop over material properties to interpolate

for vn = 1:numV
    temp = (wc1*accumarray([ICN(cell1)', ones(sum(cell1),1)], varargin{vn}(cell1).*wm1)./w1 + ...
        wc2*accumarray([ICN(cell2)', ones(sum(cell2),1)], varargin{vn}(cell2).*wm2)./w2)/...
        (wc1+wc2);
    n2interp(vn).data(2:end-1, 1) = temp(2:end);
end

%%% right edge

rightEdge = jcn==Nx-1 & icn>1 & icn<Ny;
shifted =   jcn==Nx-1 & icn<Ny-1;

% cell 1

cell1 = shifted & quad==3;

ddx = dx(Nx-1);
ddy = dy(ICN(cell1)-1);
dxm = xm(cell1) - x(Nx-1);
dym = ym(cell1) - y(ICN(cell1));
wm1 = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w1 = accumarray([ICN(cell1)', ones(sum(cell1),1)], wm1);

% cell 2

cell2 = rightEdge & quad==2;

ddx = dx(Nx-1);
ddy = dy(ICN(cell2));
dxm = xm(cell2) - x(Nx-1);
dym = ym(cell2) - y(ICN(cell2));
wm2 = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)./(ddx.*ddy);
w2 = accumarray([ICN(cell2)', ones(sum(cell2),1)], wm2);

%loop over material properties to interpolate

for vn = 1:numV
    temp = (wc1*accumarray([ICN(cell1)', ones(sum(cell1),1)], varargin{vn}(cell1).*wm1)./w1 + ...
        wc2*accumarray([ICN(cell2)', ones(sum(cell2),1)], varargin{vn}(cell2).*wm2)./w2)/...
        (wc1+wc2);
    n2interp(vn).data(2:end-1, Nx) = temp(2:end);
end

%% CORNERS

% upper left

upperLeft = jcn==1 & icn==1 & quad==1;

ddx = dx(1);
ddy = dy(1);
dxm = xm(upperLeft) - x(1);
dym = ym(upperLeft) - y(1);
wm  = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)/(ddx*ddy);
wco = sum(wm);

for vn = 1:numV
    n2interp(vn).data(1,1) = sum(varargin{vn}(upperLeft).*wm)./wco;
end

% upper right

upperRight = icn==1 & jcn==Nx-1 & quad==2;

ddx = dx(Nx-1);
ddy = dy(1);
dxm = xm(upperRight) - x(Nx-1);
dym = ym(upperRight) - y(1);
wm  = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)/(ddx*ddy);
wco = sum(wm);

for vn = 1:numV
    n2interp(vn).data(1,Nx) = sum(varargin{vn}(upperRight).*wm)./wco;
end

% lower Right

lowerRight = icn==Ny-1 & jcn==Nx-1 & quad==3;

ddx = dx(Nx-1);
ddy = dy(Ny-1);
dxm = xm(lowerRight) - x(Nx-1);
dym = ym(lowerRight) - y(Ny-1);
wm  = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)/(ddx*ddy);
wco = sum(wm);

for vn = 1:numV
    n2interp(vn).data(Ny,Nx) = sum(varargin{vn}(lowerRight).*wm)./wco;
end

% lower left

lowerLeft = icn==Ny-1 & jcn==1 & quad==4;

ddx = dx(1);
ddy = dy(Ny-1);
dxm = xm(lowerLeft) - x(1);
dym = ym(lowerLeft) - y(Ny-1);
wm  = 1 - (dxm.*dym + (ddx-dxm).*dym + (ddy-dym).*dxm)/(ddx*ddy);
wco = sum(wm);

for vn = 1:numV
    n2interp(vn).data(Ny,1) = sum(varargin{vn}(lowerLeft).*wm)./wco;
end


                                                                                                                       SiStER_interp_normal_nodes_to_markers.m                                                             0000640 0002626 0002066 00000001222 13122410535 017373  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [varm]=SiStER_interp_normal_nodes_to_markers(varN,xc,yc,xm,ym,icn,jcn)
% Transfer a variable from normal nodes to markers


%M=length(xm);
%varm=zeros(size(xm));
Nx=length(xc)+1;
Ny=length(yc)+1;
[m,~] = size(varN);

ic = icn;
jc = jcn;

jc(xm>xc(jc)) = jc(xm>xc(jc))+1;
jc(jc<2) = 2;
jc(jc>Nx-1) = Nx-1;

ic(ym>yc(ic)) = ic(ym>yc(ic))+1;
ic(ic<2) = 2;
ic(ic>Ny-1) = Ny-1;

xNnodes = [xc(jc-1); xc(jc); xc(jc); xc(jc-1)];
yNnodes = [yc(ic-1); yc(ic-1); yc(ic); yc(ic)];

IND = sub2ind(size(varN), ic, jc);
VARnodes = [varN(IND); varN(IND+m); varN(IND+1+m); varN(IND+1)];

varm = SiStER_interp_grid_to_marker_vector(xNnodes,yNnodes,VARnodes,xm,ym);

                                                                                                                                                                                                                                                                                                                                                                              SiStER_interp_normal_to_shear_nodes.m                                                               0000640 0002626 0002066 00000004411 13122410535 017034  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   %=========================================================================
function [varS]=SiStER_interp_normal_to_shear_nodes(varN,dx,dy)
% Interpolates values on normal nodes (varN) to shear nodes
% Potential shortcoming:  interpolation to side shear nodes are only
% based on the closest normal nodes and thus values will be horizontally
% uniform with 1/2 cell of the left/right sides and vertically uniform
% with 1/2 cell of the top and bottoms
%
% G.Ito 8/2016
%=========================================================================
varS=zeros(size(varN));
[Ny,Nx]=size(varS);
i1=2:Ny-1; j1=2:Nx-1; 
i2=3:Ny;   j2=3:Nx;
dydx=varS;
dydx(2:end,2:end)=dy'*dx;

%--------------------------------------------------------------------------
% interior nodes
%--------------------------------------------------------------------------
varS(i1,j1)=(varN(i1,j1).*dydx(i2,j2)/4 + varN(i1,j2).*dydx(i2,j1)/4 +...
             varN(i2,j1).*dydx(i1,j2)/4 + varN(i2,j2).*dydx(i1,j1)/4)./...
             ((dydx(i1,j1)+dydx(i1,j2)+dydx(i2,j1)+dydx(i2,j2))./4);
                   
%--------------------------------------------------------------------------
% top and bottom, excluding corners
%--------------------------------------------------------------------------
varS(1,j1)= (varN(2,j1).*dydx(2,j2)/4 + varN(2,j2).*dydx(2,j1)/4)./...
            (dydx(2,j2)/4 + dydx(2,j1)/4);
                       
varS(Ny,j1)=(varN(Ny,j1).*dydx(Ny,j2)/4 + varN(Ny,j2).*dydx(Ny,j1)/4)./...
            (dydx(Ny,j2)/4 + dydx(Ny,j1)/4);
                        
%--------------------------------------------------------------------------
% left and right, excluding corners
%--------------------------------------------------------------------------
varS(i1,1)= (varN(i1,2).*dydx(i2,2)/4 + varN(i2,2).*dydx(i1,2)/4)./...
            (dydx(i1,2)/4 + dydx(i2,2)/4);
        
varS(i1,Nx)=(varN(i1,Nx).*dydx(i2,Nx)/4 + varN(i2,Nx).*dydx(i1,Nx)/4)./...
            (dydx(i2,Nx)./4 +dydx(i1,Nx)/4);
        
%--------------------------------------------------------------------------
% corners
%--------------------------------------------------------------------------
varS(1,1) = varN(2,2); varS(1,Nx) =varN(2,Nx);
varS(Ny,1)=varN(Ny,2); varS(Ny,Nx)=varN(Ny,Nx);
                   
                   
                      
return
                                                                                                                                                                                                                                                       SiStER_interp_shear_nodes_to_markers.m                                                              0000640 0002626 0002066 00000000666 13122410535 017220  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [varm]=SiStER_interp_shear_nodes_to_markers(varS,x,y,xm,ym,icn,jcn)
% Transfer a variable from shear nodes to markers

[m, ~] = size(varS);

xnodes = [x(jcn); x(jcn+1); x(jcn+1); x(jcn)];
ynodes = [y(icn); y(icn); y(icn+1); y(icn+1)];

INDEX = sub2ind(size(varS), icn, jcn);
VARnodes = [varS(INDEX); varS(INDEX+m); varS(INDEX+m+1); varS(INDEX+1)];

varm = SiStER_interp_grid_to_marker_vector(xnodes,ynodes,VARnodes,xm,ym);

end

                                                                          SiStER_interp_shear_to_normal_nodes.m                                                               0000640 0002626 0002066 00000000563 13122410535 017040  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   %==============================================================
function [varN]=SiStER_interp_shear_to_normal_nodes(varS)
% 
% G.Ito 8/2016
%===============================================================

varN=zeros(size(varS));
varN(2:end,2:end)=(varS(1:end-1,1:end-1) + varS(2:end,1:end-1)+....
                   varS(1:end-1,2:end)   + varS(2:end,2:end))./4;
return
                                                                                                                                             SiStER_interp_velocities_to_shear_nodes_backstop.m                                                  0000640 0002626 0002066 00000004634 13122410535 021607  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   %=========================================================================
% Computes velocities at cell corners (shear nodes)
% G.Ito 7/15
%=========================================================================


%% ------------------------------------------------------------------------
% Vx is defined on the x values of cell edges so only need to interpolate
% to y values of cell tops and bottoms
% -------------------------------------------------------------------------
vxc=zeros(Ny,Nx);
% internal nodes
vxc(2:Ny-1,:)=vx(1:Ny-2,:).*(((1-dy(1:Ny-2)./(dy(1:Ny-2)+dy(2:Ny-1))))'*ones(1,Nx))+....
              vx(2:Ny-1,:).*(((1-dy(2:Ny-1)./(dy(1:Ny-2)+dy(2:Ny-1))))'*ones(1,Nx));


% Top 
if (BC.top(1)==1);  %free slip
    vxc(1,:)=vx(2,:);
elseif (BC.top(1)==0); %No slip
    vxc(1,:)=0;
else;
    disp('Problems in interp_velocities_to_shear_nodes: ');
    disp('BC.top(1) ~= 0 or 1');
    halt
end;
% Bottom
if (BC.bot(1)==1);  %free slip
    vxc(Ny,:)=vx(Ny-1,:);
    
elseif (BC.bot(1)==0); %No slip
    if (length(BC.bot)==4);
        jj=find(x < BC.bot_xbackstop(1));
        vxc(Ny,jj)=0.0;
        jj=find(x > BC.bot_xbackstop(2));
        vxc(Ny,jj)=BC.bot(4);
        jj=find(BC.bot_xbackstop(1) <= x & x <= BC.bot_xbackstop(2));
        vxc(Ny,jj)=(x(jj)-BC.bot_xbackstop(1))./(BC.bot_xbackstop(2)-BC.bot_xbackstop(1)).*BC.bot(4);
    else
      vxc(Ny,:)=0;
    end;
else;  
    disp('Problems in interp_velocities_to_shear_nodes: ');
    disp('BC.bot(1) ~= 0 or 1');
    halt
end;



%% ------------------------------------------------------------------------
% Vy is defined on y values of cell tops and bottoms so only need to 
% interpolate to x-values of cell edges
% -------------------------------------------------------------------------
vyc=zeros(Ny,Nx);
% internal values
vyc(:,2:Nx-1)=vy(:,1:Nx-2).*(ones(Ny,1)*(1-dx(1:Nx-2)./(dx(1:Nx-2)+dx(2:Nx-1))))+....
              vy(:,2:Nx-1).*(ones(Ny,1)*(1-dx(2:Nx-1)./(dx(1:Nx-2)+dx(2:Nx-1))));

%Left
if (BC.left(1)==1); %free slip
    vyc(:,1)=vy(:,2);
elseif (BC.left(1)==0); % no slip
    vyc(:,1)=0;
else;  
    disp('Problems in interp_velocities_to_shear_nodes: ');
    disp('BC.left(1) ~= 0 or 1');
    halt
end;

%Right
if (BC.right(1)==1); %free slip
    vyc(:,Nx)=vy(:,Nx-1);
    
elseif (BC.right(1)==0); % no slip
    vyc(:,Nx)=0.0;
    
else;  
    disp('Problems in interp_velocities_to_shear_nodes: ');
    disp('BC.right(1) ~= 0 or 1');
    halt
end; 
    


%

                                                                                                    SiStER_interp_velocities_to_shear_nodes.m                                                           0000640 0002626 0002066 00000004117 13122410535 017715  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   %=========================================================================
% Computes velocities at cell corners (shear nodes)
% G.Ito 7/15
%=========================================================================


%% ------------------------------------------------------------------------
% Vx is defined on the x values of cell edges so only need to interpolate
% to y values of cell tops and bottoms
% -------------------------------------------------------------------------
vxc=zeros(Ny,Nx);
% internal nodes
vxc(2:Ny-1,:)=vx(1:Ny-2,:).*(((1-dy(1:Ny-2)./(dy(1:Ny-2)+dy(2:Ny-1))))'*ones(1,Nx))+....
              vx(2:Ny-1,:).*(((1-dy(2:Ny-1)./(dy(1:Ny-2)+dy(2:Ny-1))))'*ones(1,Nx));


% Top 
if (BC.top(1)==1);  %free slip
    vxc(1,:)=vx(2,:);

elseif (BC.top(1)==0); %No slip
    vxc(1,:)=BC.top_profile;
else
    disp('Problems in interp_velocities_to_shear_nodes: ');
    disp('BC.top(1) ~= 0 or 1');
    halt
end;
% Bottom
if (BC.bot(1)==1);  %free slip
    vxc(Ny,:)=vx(Ny-1,:);
    
elseif (BC.bot(1)==0); %No slip
    vxc(Ny,:)=BC.bot_profile;

else 
    disp('Problems in interp_velocities_to_shear_nodes: ');
    disp('BC.bot(1) ~= 0 or 1');
    halt
end;



%% ------------------------------------------------------------------------
% Vy is defined on y values of cell tops and bottoms so only need to 
% interpolate to x-values of cell edges
% -------------------------------------------------------------------------
vyc=zeros(Ny,Nx);
% internal values
vyc(:,2:Nx-1)=vy(:,1:Nx-2).*(ones(Ny,1)*(1-dx(1:Nx-2)./(dx(1:Nx-2)+dx(2:Nx-1))))+....
              vy(:,2:Nx-1).*(ones(Ny,1)*(1-dx(2:Nx-1)./(dx(1:Nx-2)+dx(2:Nx-1))));

%Left
if (BC.left(1)==1); %free slip
    vyc(:,1)=vy(:,2);
elseif (BC.left(1)==0); % no slip
    vyc(:,1)=0;
else;  
    disp('Problems in interp_velocities_to_shear_nodes: ');
    disp('BC.left(1) ~= 0 or 1');
    halt
end;

%Right
if (BC.right(1)==1); %free slip
    vyc(:,Nx)=vy(:,Nx-1);
    
elseif (BC.right(1)==0); % no slip
    vyc(:,Nx)=0.0;
    
else;  
    disp('Problems in interp_velocities_to_shear_nodes: ');
    disp('BC.right(1) ~= 0 or 1');
    halt
end; 
    


%

                                                                                                                                                                                                                                                                                                                                                                                                                                                 SiStER_locate_markers_in_grid.m                                                                     0000640 0002626 0002066 00000002524 13122415305 015600  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [quad,icn,jcn] = SiStER_locate_markers_in_grid(xm,ym,x,y,dx,dy)
% [quad,icn,jcn] = SiStER_locate_markers_in_grid(xm,ym,x,y,dx,dy)
% Tells a marker which cell (and which quadrant of that cell) it belongs to.
    % icn,jcn are the indexes of the upper-left shear node of the cell
    % that a given marker is currently in
    % quad is the quadrant of the cell that contains the marker
    % (quad = 1 means bottom-right, then numbered clockwise)
    % sped up by B. Klein in Fall 2016 by using interp1 function


%% Determine Location of Markers and Quadrant of Element 
M=length(xm);
icn = zeros(1,M);
jcn = icn;
quad = icn; % quadrant 1 = bottom-right, numbered clockwise

indX = 1:length(x);
indY = 1:length(y);

Ix = interp1(x, indX, xm, 'nearest', 'extrap');
Iy = interp1(y, indY, ym, 'nearest', 'extrap');

% [~, Ix] = min(abs(bsxfun(@minus, xm, x')));
% [~, Iy] = min(abs(bsxfun(@minus, ym, y')));

jcn(xm>x(Ix))  = Ix(xm>x(Ix));
jcn(xm<=x(Ix)) = Ix(xm<=x(Ix))-1;

icn(ym>y(Iy))  = Iy(ym>y(Iy));
icn(ym<=y(Iy)) = Iy(ym<=y(Iy))-1;

jcn(jcn==0) = 1;
jcn(jcn>length(dx)) = length(dx);

icn(icn==0) = 1;
icn(icn>length(dy)) = length(dy);

disx = abs((xm-x(jcn))./dx(jcn));
disy = abs((ym-y(icn))./dy(icn));


xRIGHT = disx > 0.5;
yUP = disy > 0.5;

quad(xRIGHT & yUP) = 3;
quad(xRIGHT & ~yUP)  = 2;
quad(~xRIGHT & yUP)  = 4;
quad(~xRIGHT & ~yUP)   = 1;
                                                                                                                                                                            SiStER_MAIN_bending_beam.m                                                                          0000640 0002626 0002066 00000012235 13122417341 014312  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   % SiStER_MAIN.m
%
% Simple Stokes with Exotic Rheologies
%
% Main routine doing initialization, time loop and outputs
%
%
% J.-A. Olive, B.Z. Klein, E. Mittelstaedt, November 2014
% S.M. Howell June 2015
% G.Ito 2015-2016

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
%Put strain-rate invariant on markers so they will be advected
% might be more accurate to have two arrays of epsIIm, one interpolated
% from epsII_s and the other from epsII_m?
%----------------------------------------------------------------------------
if (PARAMS.YNPlas==1 || nonNewtonian==1)
    epsIIm=SiStER_interp_shear_nodes_to_markers(epsII_s,x,y,xm,ym,icn,jcn);
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

if exist('isee_timesteps','var')
    if (mod(t,isee_timesteps)==0);
        [n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,ep);
        aps  = n2interp(1).data./MAT(2).ecrit;
        [n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,sxym);
        SXY  = n2interp(1).data;
        [n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,sxxm);
        SXX  = n2interp(1).data;
        See_bending_beam;
        disp('Hit return to continue');
        pause;
    end
end
     

%if (mod(t,5)==0); SiStER_see; disp('Hit any key to continue'); pause; end;

end

toc;

    
                                                                                                                                                                                                                                                                                                                                                                   SiStER_MAIN.m                                                                                       0000640 0002626 0002066 00000012344 13275723525 011656  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   % SiStER_MAIN.m
%
% Simple Stokes with Exotic Rheologies
%
% Main routine doing initialization, time loop and outputs
%
%
% J.-A. Olive, B.Z. Klein, E. Mittelstaedt, November 2014
% S.M. Howell June 2015
% G.Ito 2015-2016

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
%    epsIIm=SiStER_interp_shear_nodes_to_markers(epsII_s,x,y,xm,ym,icn,jcn);
     EXYm=SiStER_interp_shear_nodes_to_markers(EXY,x,y,xm,ym,icn,jcn);
     EXXm=SiStER_interp_normal_nodes_to_markers(EXX,xc,yc,xm,ym,icn,jcn);
     if (isfield(MAT,'psi')~=0)
        dum=SiStER_interp_normal_nodes_to_markers(Rn,xc,yc,xm,ym,icn,jcn);
        epsIIm=sqrt(0.5*EXXm.^2+0.5.*(dum-EXXm).^2+EXYm.^2 + 0.25.*dum);
     else
        epsIIm=sqrt(EXXm.^2+EXYm.^2); 
     end
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
if (exist('topo_y','var'))
if (min(topo_y)<y(2))
    disp('STOPPING: Wedge Surface at top of box')
    break
end
end

end

toc;

    
                                                                                                                                                                                                                                                                                            SiStER_MAIN_ZeroDtest.m                                                                             0000640 0002626 0002066 00000012103 13122426201 013630  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   % SiStER_MAIN.m
%
% Simple Stokes with Exotic Rheologies
%
% Main routine doing initialization, time loop and outputs
%
%
% J.-A. Olive, B.Z. Klein, E. Mittelstaedt, November 2014
% S.M. Howell June 2015
% G.Ito 2015-2016

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

%--------------------------------------------------------------------------
% This is added code for Zero-D Test
%--------------------------------------------------------------------------
SXX=sxxOLD+dsxx;
SXY=sxyOLD+dsxy;
SXY_n=SiStER_interp_shear_to_normal_nodes(SXY);
SII=sqrt(SXX(2:end,2:end).^2+SXY_n(2:end,2:end).^2);
if (t==1);
    SIIout=[time mean(SII(:))];
else
    SIIout=[SIIout; [time mean(SII(:))];];
end

%----------------------------------------------------------------------------
%Put strain-rate invariant on markers so they will be advected
% might be more accurate to have two arrays of epsIIm, one interpolated
% from epsII_s and the other from epsII_m?
%----------------------------------------------------------------------------
if (PARAMS.YNPlas==1 || nonNewtonian==1)
    epsIIm=SiStER_interp_shear_nodes_to_markers(epsII_s,x,y,xm,ym,icn,jcn);
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
if (PARAMS.YNElast==1) %added or YNPlas G.Ito
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

end

toc;

    
                                                                                                                                                                                                                                                                                                                                                                                                                                                             SiStER_material_props_on_nodes.m                                                                    0000640 0002626 0002066 00000007605 13275724316 016042  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   %=========================================================================
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
    psi_s  = n2interp(1).data;  %this is needed only for epsII_s
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
end                                                                                                                           SiStER_output.m                                                                                     0000640 0002626 0002066 00000004137 13275234360 012465  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   %%
%   SiStER_output
% 
% Writes visualization output every dt_out steps as well as all data
% every dt_out_restart
%
% G.Ito 7/15
%%-----------------------------------------------------------------------

filename=num2str(t);
% if exist('dt_out_restart','var')   %This doesnt work yet because arrays
%                                    %read in here overwrite those in restart
%                                    %input file
 if (mod(t,dt_out_restart)==0 || t==Nt);
%       tout=t-floor(t/(4*dt_out_restart))*(4*dt_out_restart);
%       if tout==0; tout=4*dt_out_restart; end;
      tout=t;
      disp(['>>>>> Writing Restart output, t=' num2str(t) ' <<<<<<<<<<<'])
      save([output_dir '/restart_' num2str(tout)]);
      
 end;
% end;

if ((mod(t,dt_out)==0 && dt_out>0) || t==Nt) 

    disp(['>>>>> Writing Visualization output: ' output_dir '/' filename '  <<<<<<<<<<<'])
    [n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,ep);
    aps  = n2interp(1).data./MAT(2).ecrit;
    [n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,sxym);
    SXY  = n2interp(1).data;
    [n2interp] = SiStER_interp_markers_to_normal_nodes(xm,ym,icn,jcn,x,y,sxxm);
    SXX  = n2interp(1).data;
    if (exist('isurface','var') && isurface==1)
    save([output_dir '/' filename],'X','Y','vxc','vyc','vx','vy','p','time','dt_m','t',....
        'xm','ym','xc','yc','etam','rhom','BC','MAT','Tm','im','idm', 'PARAMS',....
        'EXX','EXX_s','EXY','EXY_n','etas','etan','ep','icn','jcn','x','y',....
        'SXX','SXY','sxxm','sxym','Zn','Zs','stratam','topo_x','phase_s','phase_n',....
        'topo_y','qd','GEOM','outit','Cohes_n','Mu_n','Cohes_s','Mu_s','ps','yield_n');
    else
    save([output_dir '/' filename],'X','Y','vxc','vyc','vx','vy','p','time','dt_m','t',....
        'xm','ym','xc','yc','etam','rhom','BC','MAT','Tm','im','idm', 'PARAMS',....
        'EXX','EXX_s','EXY','EXY_n','etas','etan','ep','icn','jcn','x','y',....
        'SXX','SXY','sxxm', 'sxym', 'Zn','Zs','stratam','qd','GEOM','outit','Cohes_n',....
        'Mu_n','Cohes_s','Mu_s','ps','yield_n','phase_s','phase_n');
    end

    
end

                                                                                                                                                                                                                                                                                                                                                                                                                                 SiStER_patch_marker_holes.m                                                                         0000640 0002626 0002066 00000013575 13274273470 014771  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [xm,ym,im,Ifix,mp,ep,idm,Tm,sxxm,sxym,pm,epsIIm,stratam]=SiStER_patch_marker_holes(icn,jcn,quad,Nx,Ny,Mquad,Mquad_crit,xm,ym,x,y,dx,dy,im,ep,idm,Tm,sxxm,sxym,pm,epsIIm,stratam)
% [xm, ym, im, Ifix, mp, ep, idm, Tm]=SiStER_patch_marker_holes(icn,jcn,quad,Nx,Ny,Mquad,Mquad_crit,xm,ym,x,y,dx,dy,im,ep,idm, Tm)
%
% seeds new markers in all quadrants where marker density has fallen below
% threshold 
% AND THEN assigns those new markers a marker index, a marker phase and a plastic
% strain based on the properties of the markers surrounding them
% (these parameters are never passed on the grid)
%
% J.-A. Olive, and B.Z. Klein, 2012-2014
%
% added Tm as input and output JAO March 30, 2015
% added stratm for passive markers of deformation, G.Ito 7/15



M=length(xm);
md_crit=Mquad_crit;
mark_per_quad=Mquad;
mdx=min(dx)/2;
mdy=min(dy)/2;

%%%%%%%%%%%%%%% LOOK FOR EMPTY (no more markers) QUADRANTS

mp = accumarray({icn, jcn, quad}, 1, [Ny-1, Nx-1, 4]);
    
% empty = find(mp==0);
% if(~isempty(empty))
%     [iEmpty,jEmpty,kEmpty] = ind2sub(size(mp), empty);
%     for n = 1:min(length(iEmpty), 100)%length(iEmpty)
%         disp('EMPTY QUADRANT: i, j, quad =')
%         disp([iEmpty(n), jEmpty(n), kEmpty(n)])
%     end
% end
    
    
%%%%% LOCATE QUADRANTS WHERE MARKER DENSITY IS BELOW THRESHOLD

     mpCrInd = find(mp<md_crit);
    [iicr, jjcr, qqcr] = ind2sub(size(mp), mpCrInd);
    iicr = iicr';
    jjcr = jjcr';
    qqcr = qqcr';
    
% NEED TO SEED NEW MARKERS IN THOSE QUADRANTS 
% SO THAT WE ARE BACK TO THE INITIAL MARKER DENSITY IN THOSE QUADRANTS
    
xrsd=[];
yrsd=[];
im_fix=[]; % marker phase
ep_fix=[]; % plastic strain
%id_fix=[]; % marker index
te_fix=[]; % temperature
sxx_fix=[]; % stress
sxy_fix=[]; % stress
epsIIm_fix=[]; %ductile part of strain invariant G.Ito
strata_fix=[]; %strata G.Ito
pm_fix=[];  %5/18
   
if ~isempty(iicr) % if there are critical quadrants
         
    for c=1:length(iicr) % go through all critical quadrants

        % the upper-left node of the corresponding cell is
        icell=iicr(c);
        jcell=jjcr(c);
        qcell=qqcr(c);
        % the quadrant area in that cell is 
        qsize=0.25*dx(jcell)*dy(icell);
        % if the smallest quadrant in the domain (area mdx*mdy) has
        % mark_per_quad markers, then this critical quadrant needs
        Nfix=ceil(qsize/(mdx*mdy))*mark_per_quad;
        Nfix=max(Nfix-mp(iicr(c),jjcr(c),qqcr(c)),1);
            
        % find the coordinates of the upper-left corner of the quadrant
        if qcell==1 || qcell==4
            xcorn=x(jcell);
        else
            xcorn=x(jcell)+dx(jcell)/2;
        end
        if qcell==1 || qcell==2
            ycorn=y(icell);
        else
            ycorn=y(icell)+dy(icell)/2;
        end
            
        % draw random marker location
        [xmrr, ymrr]=SiStER_seed_markers_uniformly(xcorn,ycorn,dx(jcell)/2,dy(icell)/2,Nfix);
            
        xrsd=[xrsd xmrr];
        yrsd=[yrsd ymrr];
            
%%%% NOW THAT THE NEW MARKERS ARE SEEDED,
%%%%% ASSIGN PARAMETERS THAT ARE NEVER STORED ON THE EULERIAN GRID
 
% the value is assigned based on the average value of the
% markers that remain in the corresponding CELL (since there's no grid value to interpolate from)
            
% CAREFUL THIS CANNOT WORK IF WE END UP WITH AN EMPTY CELL
% if that was to happen, let's just draw "im" randomly
            
% marker identity
                       
    if isempty((ep(icn==icell & jcn==jcell)))==1  

        disp('EMPTY CELL')
        disp('SOMETHING WENT VERY WRONG')        
        im_fix=1+floor(rand(1,Nfix)*max(im)); % random phase number
        ep_fix=zeros(1,Nfix);
        %id_fix=zeros(1,Nfix);
        te_fix=zeros(1,Nfix);
        sxx_fix=zeros(1,Nfix);
        sxy_fix=zeros(1,Nfix);
        press_fix=zeros(1,Nfix);  %5/18
        epsIIm_fix_buf=zeros(1,Nfix);
        strata_fix_buf=zeros(1,Nfix);
        
        
    else
        
%         xmlocal=xm(icn==icell & jcn==jcell)';
%         ymlocal=ym(icn==icell & jcn==jcell)';
%         Tmlocal=Tm(icn==icell & jcn==jcell)';
%         FT=scatteredInterpolant(xmlocal,ymlocal,Tmlocal,'nearest');
%         temp_fix=FT(xmrr,ymrr);
        
        
        % assign the average phase of the markers that are left in the cell
        phase_fix=round(mode((im(icn==icell & jcn==jcell))));
        % assign the greatest plastic strain of the markers that are left in the cell
        strain_fix=max((ep(icn==icell & jcn==jcell)));
        % assign the average temperature of the markers that are left in
        % the cell
        temp_fix=mean((Tm(icn==icell & jcn==jcell)));
        % reassign mean stress of markers left in cell
        stress_xx_fix=mean((sxxm(icn==icell & jcn==jcell)));
        stress_xy_fix=mean((sxym(icn==icell & jcn==jcell)));
        press_fix=mean((pm(icn==icell & jcn==jcell)));  %5/18
        epsIIm_fix_buf=mean((epsIIm(icn==icell & jcn==jcell))); %Gito
        strata_fix_buf=mean((stratam(icn==icell & jcn==jcell))); %Gito

    end
    
    
    
    im_fix=[im_fix phase_fix*ones(1,Nfix)];
    ep_fix=[ep_fix strain_fix*ones(1,Nfix)];
    te_fix=[te_fix temp_fix*ones(1,Nfix)];
    sxx_fix=[sxx_fix stress_xx_fix*ones(1,Nfix)];
    sxy_fix=[sxy_fix stress_xy_fix*ones(1,Nfix)];
    pm_fix=[pm_fix press_fix*ones(1,Nfix)];    %5/18
    epsIIm_fix=[epsIIm_fix epsIIm_fix_buf*ones(1,Nfix)];
    strata_fix=[strata_fix strata_fix_buf*ones(1,Nfix)]; %G.Ito

    end
        

% NOW ASSIGN PROPERTIES TO THOSE MARKERS  
Npatch=length(xrsd);
index_fix=max(idm)+1:max(idm)+Npatch;

Ifix=M+1:1:M+Npatch; % total number of markers added to fix holes in critical quadrants
xm(Ifix)=xrsd;
ym(Ifix)=yrsd;
im(Ifix)=im_fix;
ep(Ifix)=ep_fix;
idm(Ifix)=index_fix;
Tm(Ifix)=te_fix;
sxxm(Ifix)=sxx_fix;
sxym(Ifix)=sxy_fix;
pm(Ifix)=pm_fix;  %5/18
epsIIm(Ifix)=epsIIm_fix;
stratam(Ifix)=strata_fix;  %G.Ito


fprintf('\n%d%s%d%s\n', length(Ifix), ' markers added in ', length(iicr), ' cell quadrants.')
   
else
        
Ifix=0;
           
end                                                                                                                                   SiStER_plastic_seed.m                                                                               0000640 0002626 0002066 00000001630 13122410535 013546  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   %========================================================================
% SiStER_plastic_seed
% Seed weakzone by initializing ep
% S. Howell  6/15
% G. Ito
%========================================================================

if (sum(PARAMS.plastic_seed == 'random')==6)
    seeds = (rand(1,50)*MAT(2).ecrit)/2;
    a     = ceil(ym/2e3);
    b     = ceil(xm/2e3);
    ep    = seeds(a)+seeds(b);
elseif (sum(PARAMS.plastic_seed == 'box   ')==6)
%     ep((PARAMS.seed_x-PARAMS.seed_dim) < xm & xm < (PARAMS.seed_x+PARAMS.seed_dim)....
%         & (PARAMS.seed_y-PARAMS.seed_dim) <=ym & ym <= (PARAMS.seed_y+PARAMS.seed_dim))=...
      ep(((PARAMS.seed_x - xm).^2+(PARAMS.seed_y - ym).^2) <= PARAMS.seed_dim.^2)=...
          PARAMS.seed_amp*MAT(2).ecrit;
    disp(['Seeding weak box centered on x,y, dim=' num2str(PARAMS.seed_x) ',' num2str(PARAMS.seed_y) ',' num2str(PARAMS.seed_dim)]);
    
end                                                                                                        SiStER_reshape_solver_output.m                                                                      0000640 0002626 0002066 00000000631 13122410535 015550  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [p, vx, vy]=SiStER_reshape_solver_output(S,Kc,Nx,Ny)
% [p vx vy]=HIPSTER_reshape_solver_output(S,Kc,Nx,Ny)
% reshapes the output of the Stokes solver into pressure and velocity
% arrays
% B.Z. Klein, rewritten from J.-A. Olive, Spring 2013


INDp  = 1:3:length(S);
INDvx = INDp + 1;
INDvy = INDp + 2;

p  = reshape(S(INDp), Ny, Nx)*Kc;
vx = reshape(S(INDvx), Ny, Nx);
vy = reshape(S(INDvy), Ny, Nx);
                                                                                                       SiStER_rotate_stresses.m                                                                            0000640 0002626 0002066 00000000553 13122422216 014342  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   % update elastic stresses on markers following a solve (but before advection)

[ROT]=SiStER_get_rotation_rate(vx,vy,dx,dy,BC);
[om]=SiStER_interp_shear_nodes_to_markers(ROT,x,y,xm,ym,icn,jcn);

% rotate markers
alpha=om*dt_m;
sxymtemp = sxxm.*sin(2.*alpha) + sxym.*cos(2.*alpha);
sxxm = sxxm.*(cos(alpha).^2 - sin(alpha).^2)- sxym.*sin(2.*alpha);
sxym=sxymtemp;

                                                                                                                                                     SiStER_seed_markers_uniformly.m                                                                     0000640 0002626 0002066 00000001564 13122410535 015665  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [xm, ym]=SiStER_seed_markers_uniformly(x,y,dx,dy,N)

%
% seeds N markers in the cell or quadrant whose upper-left node 
% has coordinates (x,y) and width, height (dx,dy)


% create a subgrid

fact=0.4; % randomization factor

nx=ceil(sqrt(N*dx/dy));
ny=ceil(N/nx);
actual_N=nx*ny;

ddx=linspace(x,x+dx,nx+1);
ddy=linspace(y,y+dy,ny+1);
[DX DY]=meshgrid(ddx,ddy);

xm=zeros(1,actual_N);
ym=xm;

k=0;
for i=1:ny
    for j=1:nx
        
        k=k+1;
        xsub=DX(i,j);
        ysub=DY(i,j);
        dxsub=dx/nx;
        dysub=dy/ny;
       
        % randomize marker position around sub-cell center
        xm(k)=xsub+(dxsub/2)*(1+fact*2*(rand-0.5));
        ym(k)=ysub+(dysub/2)*(1+fact*2*(rand-0.5));

    end
end


if actual_N>N % only keep N random markers out of actual_N markers

    idx=randperm(actual_N);
    idx=idx(1:N);
    xm=xm(idx);
    ym=ym(idx);
    
end                                                                                                                                            SiStER_set_decollement_props.m                                                                      0000640 0002626 0002066 00000001201 13123541765 015505  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   %=========================================================================
% Decollment is now fixed in depth.  Set properties based on:
% decollement_ytop
% decollement_mu
% decollement_coh
% Deals with material properties on grid so decollement_ytop had better
% correspond to a shear node depth
% G.Ito 4/5/17
%==========================================================================
Mu_s(y>=PARAMS.decollement_ytop,:)=PARAMS.decollement_mu;
Mu_n([0 yc]>=PARAMS.decollement_ytop,:)=PARAMS.decollement_mu;

Cohes_s(y>=PARAMS.decollement_ytop,:)=PARAMS.decollement_coh;
Cohes_n([0 yc]>=PARAMS.decollement_ytop,:)=PARAMS.decollement_coh;

                                                                                                                                                                                                                                                                                                                                                                                               SiStER_set_timestep.m                                                                               0000640 0002626 0002066 00000000620 13122410535 013612  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [dt_m]=SiStER_set_timestep(dx,dy,vx,vy,PARAMS)
% [dt_m]=SiStER_set_timestep(dx,dy,vx,vy,PARAMS)
% sets the advection time step
% J.-A. Olive, November 2014

%dt_m=PARAMS.fracCFL*0.5*min(min(dx),min(dy))./max(max(max(abs(vx))),max(max(abs(vy))));
dt1=min(dy'./max(abs(vy(1:end-1,:)'+vy(2:end,:)')/2)');
dt2=min(dx./max(abs(vx(:,2:end)+vx(:,1:end-1))./2));
dt_m=PARAMS.fracCFL.*min([dt1 dt2]);                                                                                                                SiStER_thermal_solver_sparse.m                                                                      0000640 0002626 0002066 00000013657 13122410535 015526  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   function [T, rhs, Lii, Ljj, Lvv]=SiStER_thermal_solver_sparse(x,y,Told,rho,cp,kx,ky,dt,BCtherm)
% implicit Solver for thermal diffusion
% rho cp dT/dt = div (k grad T)
% J.-A. Olive, November 2014
% B.Z. Klein, added sparse matrix filling, End of 2014


% need to use sparse-filling for L

Nx=length(x);
Ny=length(y);
dx=diff(x);
dy=diff(y);

Li=zeros(Nx*Ny, 1);
Lj=zeros(Nx*Ny, 1);
Lv=zeros(Nx*Ny, 1);
rhs=zeros(Nx*Ny,1);

n = 1;

for i=1:Ny
    for j=1:Nx
        
        in=i+Ny*(j-1); % global index
        
        if i==1 % top boundary
            
            if BCtherm.top(1)==1 % Dirichlet
            
                Li(n) = in;
                Lj(n) = in;
                Lv(n)  = 1;
                
                n = n+1;
               
                rhs(in) = BCtherm.top(2);
                
            %L(in,in)=1;
            %rhs(in)=BCtherm.top(2);
            
            elseif BCtherm.top(1)==2 % Neumann
            
                Li(n) = in;
                Lj(n) = in;
                Lv(n) = -1;
                
                n = n+1;
                
                Li(n) = in;
                Lj(n) = in+1;
                Lv(n) = 1;
                
                n = n+1;
                
                rhs(in)=BCtherm.top(2)*dy(i);
                
%             L(in,in)=-1;
%             L(in,in+1)=1;
%             rhs(in)=BCtherm.top(2)*dy(i);    
            
            end
            
        elseif i==Ny % bottom boundary
            
            
            if BCtherm.bot(1)==1 % Dirichlet
            
                Li(n) = in;
                Lj(n) = in;
                Lv(n) = 1;
                
                n = n+1;
                
                rhs(in) = BCtherm.bot(2);
                
%             L(in,in)=1;
%             rhs(in)=BCtherm.bot(2);
            
            elseif BCtherm.bot(1)==2 % Neumann
                
                Li(n) = in;
                Lj(n) = in;
                Lv(n) = 1;
                
                n = n+1;
                
                Li(n) = in;
                Lj(n) = in-1;
                Lv(n) = -1;
                
                n = n+1;
                
                rhs(in) = BCtherm.bot(2)*dy(i-1);
                
%             L(in,in)=1;
%             L(in,in-1)=-1;
%             rhs(in)=BCtherm.bot(2)*dy(i-1);    
            
            end
 
        elseif j==1 % left boundary
            
            
            if BCtherm.left(1)==2 % Neumann
            
                Li(n) = in;
                Lj(n) = in;
                Lv(n) = -1;
                
                n = n+1;
                
                Li(n) = in;
                Lj(n) = in+Ny;
                Lv(n) = 1;
                
                n = n+1;
            
                rhs(in) = BCtherm.left(2)*dx(j);
                
%             L(in,in)=-1;
%             L(in,in+Ny)=1;
%             rhs(in)=BCtherm.left(2)*dx(j);  
            
            elseif BCtherm.left(1)==1 % Dirichlet
            
                Li(n) = in;
                Lj(n) = in;
                Lv(n) = 1;
                
                n = n+1;
                
                rhs(in) = BCtherm.left(2);
                
%             L(in,in)=1;
%             rhs(in)=BCtherm.left(2); 
            
            end
            
            
        elseif j==Nx % right boundary
            
            
            if BCtherm.right(1)==2 % Neumann
            
                Li(n) = in;
                Lj(n) = in;
                Lv(n) = 1;
                
                n = n+1;
                
                Li(n) = in;
                Lj(n) = in-Ny;
                Lv(n) = -1;
                
                n = n+1;
                
                rhs(in) = BCtherm.right(2)*dx(j-1);
                
%             L(in,in)=1;
%             L(in,in-Ny)=-1;
%             rhs(in)=BCtherm.right(2)*dx(j-1);  
            
            elseif BCtherm.right(1)==1 % Dirichlet
                
                Li(n) = in;
                Lj(n) = in;
                Lv(n) = 1;
                
                n = n+1;
                
                rhs(in) = BCtherm.right(2);
                
%             L(in,in)=1;
%             rhs(in)=BCtherm.right(2); 
            
            end
            
            
            
        else
            
            
%         %   internal nodes
            
            ddx=dx(j-1)+dx(j);
            ddy=dy(i-1)+dy(i);
            
            Li(n) = in;
            Lj(n) = in;
            Lv(n) = rho(i,j)*cp(i,j)+2*dt*kx(i,j)/(dx(j)*ddx) + 2*dt*kx(i,j-1)/(dx(j-1)*ddx) + ...
                2*dt*ky(i,j)/(ddy*dy(i)) + 2*dt*ky(i-1,j)/(ddy*dy(i-1));
            
            n = n+1;
            
%             L(in,in)=rho(i,j)*cp(i,j)+2*dt*kx(i,j)/(dx(j)*ddx) + 2*dt*kx(i,j-1)/(dx(j-1)*ddx) + ...
%                 2*dt*ky(i,j)/(ddy*dy(i)) + 2*dt*ky(i-1,j)/(ddy*dy(i-1));
            
            Li(n) = in;
            Lj(n) = in+Ny;
            Lv(n) = -2*dt*kx(i,j)/(dx(j)*ddx);
            
            n = n+1;
            
%             L(in,in+Ny)=-2*dt*kx(i,j)/(dx(j)*ddx);
% 
            Li(n) = in;
            Lj(n) = in-Ny;
            Lv(n) = -2*dt*kx(i,j-1)/(dx(j-1)*ddx);
            
            n = n+1;
% 
% %             L(in,in-Ny)=-2*dt*kx(i,j-1)/(dx(j-1)*ddx);
%             
            Li(n) = in;
            Lj(n) = in+1;
            Lv(n) = -2*dt*ky(i,j)/(dy(i)*ddy);
            
            n = n+1;
% 
% %             L(in,in+1)=-2*dt*ky(i,j)/(dy(i)*ddy);
% 
            Li(n) = in;
            Lj(n) = in-1;
            Lv(n) = -2*dt*ky(i-1,j)/(dy(i-1)*ddy);
            
            n = n+1;

%             L(in,in-1)=-2*dt*ky(i-1,j)/(dy(i-1)*ddy);
            
            rhs(in)=rho(i,j)*cp(i,j)*Told(i,j);
        
        end
        
    end
end

nn = n-1;

Lii = Li(1:nn);
Ljj = Lj(1:nn);
Lvv = Lv(1:nn);

L = sparse(Lii, Ljj, Lvv);

tvec=L\rhs;

T=reshape(tvec,Ny,Nx);

    

    
        
        
        
        
                                                                                 SiStER_thermal_update.m                                                                             0000640 0002626 0002066 00000002750 13122410535 014111  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   % SiStER THERMAL SOLVE

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
Tm=Tm+dTm;                        SiStER_update_ep.m                                                                                  0000640 0002626 0002066 00000003324 13275724524 013076  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   
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
% epsII_s (used in the old method above) in flow_solve.  These are strain changes, 
% but using (unneeded) arrays for strainrate to save RAM

EXY=0.*EXY; EXX=EXY;  Rn=0.*Rn;
EXY(s_nodes_yield)=EXY(s_nodes_yield).*dt_m;
EXX(n_nodes_yield)=EXX(n_nodes_yield).*dt_m;

EXYm=SiStER_interp_shear_nodes_to_markers(EXY,x,y,xm,ym,icn,jcn);  
EXXm=SiStER_interp_normal_nodes_to_markers(EXX,xc,yc,xm,ym,icn,jcn);

if (isfield(MAT,'psi'));
    Rn=Rn.*dt_m;
    dum=SiStER_interp_normal_nodes_to_markers(Rn,xc,yc,xm,ym,icn,jcn);
    dum=sqrt(0.5.*EXXm.^2+0.5*(dum-EXXm).^2+EXYm.^2 + 0.25.*dum);
else
    dum=sqrt(EXXm.^2+EXYm.^2);
end


ep=(ep+dum)./(dt_m/PARAMS.tau_heal+1);
                                                                                                                                                                                                                                                                                                            SiStER_update_markers.m                                                                             0000640 0002626 0002066 00000005003 13274273563 014133  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   % SiStER Update Markers
% Added stratam for passively visualizing deformation G.Ito 7/18
% Added pm G.Ito 5/18


[xm_new,ym_new] = SiStER_advect_markers(x,y,xm,ym,dx,dy,dt_m,vx,vy);
xm=xm_new;
ym=ym_new;
% eliminate markers that left domain
Iin=find(xm<=xsize & xm>=0 & ym>=0 & ym<=ysize);

msg2='  markers removed: ';
msg=[msg2 num2str(length(xm)-length(Iin))];
disp(msg)
xm=xm(Iin);
ym=ym(Iin);
im=im(Iin);
ep=ep(Iin);
Tm=Tm(Iin);
idm=idm(Iin);
sxxm=sxxm(Iin);
sxym=sxym(Iin);
epsIIm=epsIIm(Iin);
stratam=stratam(Iin);
pm=pm(Iin);  %5/18


% locate advected markers with respect to the eulerian grid
[quad,icn,jcn] = SiStER_locate_markers_in_grid(xm,ym,x,y,dx,dy);
    
% check for holes in the marker distribution, 
% patch with new markers if necessary
% those new markers immediately get assigned a value of phase (im), index 
% (idm) and accumulated plastic strain (ep), i.e., the 2 variables that never get
% passed to nodes. 
[xm, ym, im, Ifix, mp, ep, idm, Tm, sxxm, sxym, pm, epsIIm, stratam]=...
    SiStER_patch_marker_holes(icn,jcn,quad,Nx,Ny,Mquad,Mquad_crit,xm,ym,x,y,dx,dy,im,ep,idm,Tm,sxxm,sxym,pm,epsIIm,stratam);   %G.Ito

% then they get assigned P, epsII and stresses from grid values

if min(Ifix)>0
    
    msg2='  markers added: ';
    msg=[msg2 num2str(length(Ifix))];
    disp(msg)
    xmFIX=xm(Ifix);
    ymFIX=ym(Ifix);
    
    % pass temperature, pressure, strain rate and stresses to the new
    % markers from their nodal values
    
    % locate new markers with respect to the eulerian grid
    [quadFIX,icnFIX,jcnFIX] = SiStER_locate_markers_in_grid(xmFIX,ymFIX,x,y,dx,dy);
    
    
    % interpolate grid values to the new markers
    % from shear nodes:

    
     % dont do that- you'd be passing pre-advection nodal Ts to markers that were just advected    
     %[temp]=SiStER_interp_shear_nodes_to_markers(T,x,y,xmFIX,ymFIX,icnFIX,jcnFIX);
     %Tm(Ifix)=temp; % temperature
     
% UNCOMMENT WHEN I ADD  STRESSES    
%     [temp]=SiStER_interp_shear_nodes_to_markers(sxy0,x,y,xmFIX,ymFIX,icnFIX,jcnFIX);
%     sxym(Ifix)=temp; % stress_xy
%     
%     % from normal nodes
% 
%     [temp]=SiStER_interp_normal_nodes_to_markers(sxx0,xc,yc,xmFIX,ymFIX,icnFIX,jcnFIX);
%     sxxm(Ifix)=temp; % stress_xx
% removed G.Ito 5/4/18
%     [temp]=SiStER_interp_normal_nodes_to_markers(p,xc,yc,xmFIX,ymFIX,icnFIX,jcnFIX);
%     pm(Ifix)=temp; % pressure
    
    % DO IT ALSO FOR STRAIN RATE !
    
end
    
    
% locate all markers with respect to the eulerian grid

[qd,icn,jcn] = SiStER_locate_markers_in_grid(xm,ym,x,y,dx,dy);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             SiStER_update_marker_stresses.m                                                                     0000640 0002626 0002066 00000001652 13275720302 015676  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   %=============================================================================
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
pm=pm+dum;                                                                                      SiStER_update_surface.m                                                                             0000640 0002626 0002066 00000006007 13217147674 014124  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   %% ========================================================================
% SiStER_update_surface
% advects surface, erodes it, and then ensures they remain across the
% surface for contractional cases
% G. Ito modified from J. Weiss's code 11/15
% G.Ito overhanging topo. now removed BEFORE erosion calculation 7/4/17
%% ========================================================================

if (exist('isurface','var'))
if (isurface==1)
    if (BC.right(3)>=0);
        disp('Surface topography works only for inflow on left side');
        disp('Stopped in SiStER_update_surface.m')
        halt
    end
    % UPDATE MARKER CHAIN FOR TOPOGRAPHY
    [topo_x,topo_y] = SiStER_advect_markers(x,y,topo_x,topo_y,dx,dy,dt_m,vx,vy);
    %Remove any cliffs or overhanging topography
    ii=find(topo_x(2:end)-topo_x(1:end-1)<(xsize./Ntopo/1e3));
    if (length(ii)>=1);
       topo_x(ii+1)=[];
       topo_y(ii+1)=[];
    end;
    %topo_x must remain in order of increasing x. 
    [topo_x,isort]=sort(topo_x);
    topo_y=topo_y(isort);
% MAKE SURE TOPO MARKERS SPAN FULL WIDTH OF DOMAIN
% Original code for extension case
%    Iin=find(topo_x<=xsize & topo_x>=0);
%    topo_x=topo_x(min(Iin)-1:max(Iin)+1);
%    topo_y=topo_y(min(Iin)-1:max(Iin)+1);
    num=round((xsize-topo_x(end))/PARAMS.dxsurf); % number of markers to add (~10 m increment)
    if (num>0)
        new_x=linspace(topo_x(end)+PARAMS.dxsurf,xsize,num); % new x values
        new_y=topo_y(end)*ones(size(new_x)); % new y values
        topo_x=[topo_x new_x];
        topo_y=[topo_y new_y];
    end;

% Reseed topography if needed
    reseed_yes=0;
    if (PARAMS.kappa>0);
        dt_surf=0.5*min(diff(topo_x)).^2/PARAMS.kappa;
        nsolve=ceil(dt_m/dt_surf);
        if (nsolve>1e3);
            reseed_yes=1;
        end;
    end;
    if (mod(t,PARAMS.reseedsurface)==0 || reseed_yes)
        disp(['Reseeding topography t = ' num2str(t)])
        sorter = sortrows([topo_x' topo_y'],1);    % makes sure things are in the ride order
        topo_x = sorter(:,1)'; 
        topo_y = sorter(:,2)';
        %Ntopo  = 1e5;   % number of topo markers
        Ntopo=round(xsize/PARAMS.dxsurf)+1;
        dum1  = topo_x;
        topo_x =linspace(0,xsize,Ntopo);
        topo_y = interp1(dum1,topo_y,topo_x);
    end
% EROSION
tcpu4=cputime;
    if PARAMS.kappa>0    
        topo_y=erosion_by_diffusion(topo_x,topo_y,dt_surf,dt_m,PARAMS.kappa);
    
    % % RESET TOPOGRAPHY
    % % i.e. make everything above/below surface air/rock
    % % added by GI/JW/JAO - 07/15
        %ind=find(isnan(topo_y));topo_x(ind)=[];topo_y(ind)=[]; % interp1 choking on NaN/Inf values for some diffusivities - JW 07/15
        ii=find(topo_x>-0.1*xsize & topo_x<1.1*xsize & topo_y>-0.1*ysize & topo_y<1.1*ysize);
        topo_x=topo_x(ii); topo_y=topo_y(ii);
      %Removing overhanging topomarkers WAS here
        topomarkers=interp1(topo_x,topo_y,xm,'spline','extrap');
        im(ym<topomarkers & im==2)=1;
        ii=ym>=topomarkers & im==1;
        im(ii)=2;
        stratam(ii)=-2;
        
    end
end
end
tcpu5=cputime;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         SiStER_VEP_rheology.m                                                                               0000640 0002626 0002066 00000006462 13275724470 013500  0                                                                                                    ustar   gito                            gito                                                                                                                                                                                                                   %==========================================================================
% SiStER_UPDATE_RHEOLOGY
% calculates visco-elasto-plastic rheology terms on shear and normal nodes
% G.Ito 8/20
%==========================================================================

%--------------------------------------------------------------------------
% plastic viscosity first (Mohr-Coulomb)
%--------------------------------------------------------------------------
if (pit==1);  
    ps=ps_old;  %ps_old is from pm, just advected from previous t-step G.Ito 5/18
else
    dp=p-pold;  %update ps from previous Picard iteration
    ps=ps_old+SiStER_interp_normal_to_shear_nodes(dp,dx,dy);  %<<<< Interpolation done here
end

if (PARAMS.YNPlas==1);
    ps(ps<0)=0;
    pnn=p; 
    pnn(p<0)=0;
    yield_s=(Cohes_s+Mu_s.*ps).*cos(atan(Mu_s));
    yield_n=(Cohes_n+Mu_n.*pnn).*cos(atan(Mu_n));
    if (PARAMS.welast==0);
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



                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              