%=========================================================================
% Decollment is now fixed in depth.  Set properties based on:
% decollement_ytop
% decollement_mu
% decollement_coh
% Deals with material properties on grid so decollement_ytop had better
% correspond to a shear node depth
% G.Ito 4/5/17
%==========================================================================
if (isfield(PARAMS,'decollement_ratio')==0)
    Mu_s(y>=PARAMS.decollement_ytop,:)=PARAMS.decollement_mu;
    Mu_n([0 yc]>=PARAMS.decollement_ytop,:)=PARAMS.decollement_mu;
else
    Mu_s(y>=PARAMS.decollement_ytop,:)=Mu_s(y>=PARAMS.decollement_ytop,:).*PARAMS.decollement_ratio;
    Mu_n([0 yc]>=PARAMS.decollement_ytop,:)=Mu_n([0 yc]>=PARAMS.decollement_ytop,:).*PARAMS.decollement_ratio;
end

Cohes_s(y>=PARAMS.decollement_ytop,:)=PARAMS.decollement_coh;
Cohes_n([0 yc]>=PARAMS.decollement_ytop,:)=PARAMS.decollement_coh;

