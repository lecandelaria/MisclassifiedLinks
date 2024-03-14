function [V_hat] = V100_Avarmis(V110_Avarnetmis, vils, rho)
% Outputs an estimate of the asymptotic variance 
%
% vils - sample of networks.
% lambda - smoothing parameter.
% rho - vector of misclassified probabilities.

load(['results_mis/theta_hat_l0_r0_',num2str(rho(1)),'_r1_',num2str(rho(2)) ,'.mat'], 'theta_hat');

r0=rho(1);
r1=rho(2);
p = size(theta_hat,1); 
V = zeros(p,p);
Psi = zeros(p,p);

for w=vils
    
    fprintf('        village %d\n',w);
    g = csvread(['directed_adjacency_matrices/lendmoney',num2str(w),'.csv']);
    N =  size(g,1) - 1 ;
    load(['dMU/exstat_vil',num2str(w),'.mat']);
    load(['dMU/enstat_vil',num2str(w),'_lam0.mat']);
    load(['dMU/psigamma_vil',num2str(w),'_lam0.mat']);
    
    c = - cat(2, eye(3), zeros(3,1))*([1- r0-r1 , 0 , 0, 0 ; 0, 1-r0-r1, 0, 0 ; 0 , 0, (1-r0-r1)^2, r0*(1-r0-r1); 0 , 0, 0 , 1-r0-r1 ]\[r0;r0;r0^2;r0] );
    C =  cat(2, eye(3), zeros(3,1))*inv([1- r0-r1 , 0 , 0, 0 ; 0, 1-r0-r1, 0, 0 ; 0 , 0, (1-r0-r1)^2, r0*(1-r0-r1); 0 , 0, 0 , 1-r0-r1 ]);

% Computing Endogenous Statistics.
    stat_en_sum  =  stat_en(:,:,3)+ stat_en(:,:,3)';
    stat_en      =  cat(3, stat_en(:,:,1), stat_en(:,:,3), stat_en(:,:,2), stat_en_sum);
    stat_en      =  repmat( reshape(c,[1,1,3]), [ size(g,1)-1,size(g,1)-1, 1 ] ) +  multiprod(stat_en, C', [2,3], [1,2]) ;

    [V_new, Psi_new] = V110_Avarnetmis(g, stat_en, stat_ex, theta_hat, psi_gamma, rho);
    
    V = V + V_new;
    Psi = Psi + (Psi_new./N);

end

V_hat = V\Psi/V;

save(['results_mis/Avar_hat_l0_r0_',num2str(rho(1)),'_r1_',num2str(rho(2)) ,'.mat'], 'V_hat');

    
end
% Last update: Jan 2022.

