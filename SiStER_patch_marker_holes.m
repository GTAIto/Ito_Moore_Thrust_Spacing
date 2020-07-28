function [xm,ym,im,Ifix,mp,ep,idm,Tm,sxxm,sxym,pm,epsIIm,stratam]=SiStER_patch_marker_holes(icn,jcn,quad,Nx,Ny,Mquad,Mquad_crit,xm,ym,x,y,dx,dy,im,ep,idm,Tm,sxxm,sxym,pm,epsIIm,stratam)
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
           
end