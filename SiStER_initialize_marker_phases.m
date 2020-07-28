function [im] = SiStER_initialize_marker_phases(Nphase,GEOM,xm,ym,PARAMS,xsize)

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

end