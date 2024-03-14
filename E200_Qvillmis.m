function [Ln, dQn] = E200_Qvillmis(E210_dQmis, vils, rho, theta)

p = size(theta,1);
Ln  = zeros(1) ;
dQn = zeros(p,1) ;

% Sum over villages
for w=vils
    
    g = csvread(['directed_adjacency_matrices/lendmoney',num2str(w),'.csv']);
    
    load(['dMU/exstat_vil',num2str(w),'.mat']); 
    load(['dMU/enstat_vil',num2str(w),'_lam0.mat']);
      
    [L0, dQ0] = E210_dQmis(g, stat_en, stat_ex, rho, theta);
    Ln  = Ln + L0 ;
    dQn = dQn + dQ0 ; 


end
    