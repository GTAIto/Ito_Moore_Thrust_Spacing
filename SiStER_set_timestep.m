function [dt_m]=SiStER_set_timestep(dx,dy,vx,vy,PARAMS)
% [dt_m]=SiStER_set_timestep(dx,dy,vx,vy,PARAMS)
% sets the advection time step
% J.-A. Olive, November 2014

%dt_m=PARAMS.fracCFL*0.5*min(min(dx),min(dy))./max(max(max(abs(vx))),max(max(abs(vy))));
dt1=min(dy'./max(abs(vy(1:end-1,:)'+vy(2:end,:)')/2)');
dt2=min(dx./max(abs(vx(:,2:end)+vx(:,1:end-1))./2));
dt_m=PARAMS.fracCFL.*min([dt1 dt2]);