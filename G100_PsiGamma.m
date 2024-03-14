function G100_PsiGamma(vils,lambda,slow)
% Outputs and saves an [N x N x p] matrix for each village, where the ijth
% element is the vector \hat E[S_{ij} | X, \sigma]. Thus, this is the
% matrix of of derivatives of the marginal utilities V_{ij}, where the
% derivative is taken with respect to the structural parameter \theta.
%
% 'slow' is an indicator for whether or not to use en_stat_slow in place of
% en_stat

load('X.mat');

for w=vils
        
    fprintf('Psi_gamma, village=%d\n',w);
    
    g = csvread(['directed_adjacency_matrices/lendmoney',num2str(w),'.csv']);
    Z = csvread(['directed_adjacency_matrices/rel',num2str(w),'.csv']);
           
    tstart=tic;
    
    if(slow ~= 1)
        psi_gamma = G110_LinearEnstat(g, Xest, Z, lambda) ;
    else
        psi_gamma = G111_LinearEnstatSlow(g, Xest, Z, lambda) ;
    end

    telapsed = toc(tstart);

    fprintf('    Time elapsed: %.3g.\n', telapsed);
    save(['dMU/psigamma_vil',num2str(w),'_lam0.mat'], 'psi_gamma');
    

end
