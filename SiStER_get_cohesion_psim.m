function [cohes,psim]=SiStER_marker_cohesion_psim(im,ep,MAT)
% [cohes]=SiStER_get_cohesion(im,ep,MAT)
% compute cohesion on markers based on ep
%G.Ito 8/2016
% Now also computes psim for dilation
% sped up by B. Klein 9/2016
%Added ep_startweakening 4/17

cohes = zeros(size(im));
if (isfield(MAT,'psi')~=0)
    psim=zeros(size(im));
end;

types = unique(im);
for i = 1:length(types)
    imInd = find(im == types(i));
    Cmax=MAT(types(i)).Cmax;
    Cmin=MAT(types(i)).Cmin;
    epscrit=MAT(types(i)).ecrit;
    ep1=MAT(types(i)).ep_startweakening; 
    %ep(imInd)=min(ep(imInd),epscrit);
    %disp(['epscrit=',num2str(epscrit)]);
    % get cohesion
    cohes(imInd)=min( max( Cmax+(Cmin-Cmax).*(ep(imInd)-ep1)./(epscrit-ep1),Cmin),Cmax);
    if (isfield(MAT,'psi'))
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


