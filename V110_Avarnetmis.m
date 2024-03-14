function [V,Psi] = V110_Avarnetmis(g, stat_en, stat_ex, theta_hat, psi_gamma, rho)
% Outputs an estimate of the asymptotic variance. 
%
% g - adjacency matrix.
% stat_en - vector of endogenous statistics (the output of en_stat.m). 
% stat_ex - vector of exogenous statistics (the output of ex_stat.m). 
% theta_hat - structural parameter. 
% Xest - from X.mat, and 
% rho - misspecified probabilities. 
%
% This function uses multiprod.m.

G   = g(2:size(g,1),2:size(g,1));
N   = size(G,1);     % number of agents
p   = length(theta_hat);
r0  = rho(1);
r1  = rho(2);

C =  cat(2, eye(3), zeros(3,1))*inv([1- r0-r1 , 0 , 0, 0 ; 0, 1-r0-r1, 0, 0 ; 0 , 0, (1-r0-r1)^2, r0*(1-r0-r1); 0 , 0, 0 , 1-r0-r1 ]);

ES_ij = cat(3, stat_en, stat_ex);                   % [N x N x p] ijth entry = \hat E[S_{ij} | X, sigma]
V_ij  = multiprod(ES_ij, theta_hat, [2 3], [1 2]);  % [N x N] ijth entry = \hat E[S_{ij} | X, \sigma]' \theta

% Quasi-Likelihood
Phi         =  normcdf(V_ij) ; 
Phi_mis     =  r0 + (1-r0-r1).*Phi ; 
Phi_mis(Phi_mis==1)=0.999 ;
Phi_mis(Phi_mis==0)=0.001 ;
phi         =  normpdf(V_ij);
phi(1:(size(phi,1)+1):end) = 0; % zero out diagonal elements

den = Phi_mis.*(1-Phi_mis);
num = multiprod(ES_ij, (1-r0-r1).*phi , 3) ;

%% V_n (second derivative of the log-likelihood)

H1 =  (1-r0-r1).*V_ij.*phi.*(Phi_mis -G) - ((1-r0-r1).*phi).^2 ;
H2 =  (G - Phi_mis).*( ((1-r0-r1).*phi).^2 ).*(1-2.*Phi_mis) ;

H = H1./(den) - H2./(den.^2) ;

ESESt = multiprod(ES_ij, permute(ES_ij, [1 2 4 3]), [3 4]); % [N x N x p x p] equals E[S | X, \sigma] * E[S | X, \sigma]'

V_sum  = multiprod(ESESt, H , [3 4], 3);              % [N x N x p x p]
V      = squeeze(sum(sum(V_sum,1),2))/(N^2);      % [p x p]
clear ESESt V_sum

%% Psi_n = Psi_A + Psi_B
    
    Psi_A =  permute( sum(  multiprod(ES_ij, (( G - Phi_mis ).*(1-r0-r1).*phi)./den , 3) , 2)./N, [1 3 2] );
             % [N x p], the rows are k specific.
   
    % Psi_B
    %----------------------------------------------------------------------
    Dg_Psi =  reshape( ((1-r0-r1).*phi) , [N^2, 1]).*(theta_hat(1:3)'*C)      ;  % size(Dg_Psi)[N^2, 4]
    v1_ij  =  reshape( ((1-r0-r1).*phi), [N^2, 1]).*reshape( ES_ij, [N^2, p]) ;  % size(v1_ij) [N^2, p]
    
   Psi_B = zeros(N,p);
   i =1; 
   
   while i < (N^2)
       
       Psi_B_temp   = -v1_ij(i, :)'*Dg_Psi(i,:) ; 
       psigamma_temp = squeeze(psi_gamma(i, :, :) )' ;
       Psi_B = Psi_B + (Psi_B_temp*psigamma_temp)'./N^2 ;
       
       i = i +1 ;
   end
   
   Psi_k = Psi_A + Psi_B; 
   
   Psi_k = Psi_k - mean(Psi_k) ; 
   Psi = (Psi_k'*Psi_k)/N ;
   
end
       
%Last Update: Dec 20, 2019.