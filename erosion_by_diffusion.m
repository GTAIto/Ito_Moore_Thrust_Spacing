function [topo_new]=erosion_by_diffusion(xc,topo_old,dt_surf,dt,K)

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
