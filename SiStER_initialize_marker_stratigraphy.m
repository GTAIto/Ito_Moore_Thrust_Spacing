function [stratam] = SiStER_initialize_marker_strata(strata_y,xm,ym)

% assign passive stratigraphic horizons
% these do not change material properties, just used to show deformation
% G.Ito 7/15
%-------------------------------------------------------------------------
stratam=zeros(size(xm));


for kk = 1:length(strata_y)-1
    stratam(ym>=strata_y(kk) & ym< strata_y(kk+1) )=kk;
end

stratam(ym>=strata_y(end)) = length(strata_y);

end