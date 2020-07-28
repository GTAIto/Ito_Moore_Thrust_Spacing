function [xm, ym] = SiStER_initialize_marker_positions(xsize,ysize,dx,dy,Mquad)
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


