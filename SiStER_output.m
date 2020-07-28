%%
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
    
    if (exist('topo_x','var')==0); topo_x=0; topo_y=0; end;
    if (exist('psi_n','var')==0); psi_n=0; end;   
    if (exist('Plithn','var')==0); Plithn=0; Pliths=0; end;
    save([output_dir '/' filename],'PARAMS','GEOM','MAT','BC','outit',....
        'x','y','xc','yc','X','Y','vxc','vyc','vx','vy','p','time','dt_m','t',....
        'xm','ym','Tm','im','icn','jcn','qd','stratam','ep','sxxm','sxym',....
        'EXX','EXY','etas','etan','Zn','phase_s','phase_n','Cohes_s','Mu_s','yield_s',...
        'rho','Rn','psi_n','epsII_s','topo_x','topo_y','Plithn','Pliths');  
end

