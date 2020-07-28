%========================================================================
% SiStER_plastic_seed
% Seed weakzone by initializing ep
% S. Howell  6/15
% G. Ito
%========================================================================

if (sum(PARAMS.plastic_seed == 'random')==6)
    seeds = (rand(1,50)*MAT(2).ecrit)/2;
    a     = ceil(ym/2e3);
    b     = ceil(xm/2e3);
    ep    = seeds(a)+seeds(b);
elseif (sum(PARAMS.plastic_seed == 'box   ')==6)
%     ep((PARAMS.seed_x-PARAMS.seed_dim) < xm & xm < (PARAMS.seed_x+PARAMS.seed_dim)....
%         & (PARAMS.seed_y-PARAMS.seed_dim) <=ym & ym <= (PARAMS.seed_y+PARAMS.seed_dim))=...
      ep(((PARAMS.seed_x - xm).^2+(PARAMS.seed_y - ym).^2) <= PARAMS.seed_dim.^2)=...
          PARAMS.seed_amp*MAT(2).ecrit;
    disp(['Seeding weak box centered on x,y, dim=' num2str(PARAMS.seed_x) ',' num2str(PARAMS.seed_y) ',' num2str(PARAMS.seed_dim)]);
    
end