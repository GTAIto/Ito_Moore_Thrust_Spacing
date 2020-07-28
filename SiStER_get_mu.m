function [mum]=SiStER_get_mu(im,ep,MAT)
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


