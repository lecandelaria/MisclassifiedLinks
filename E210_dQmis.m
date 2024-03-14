function [L,dQ] = E210_dQmis(g, stat_en, stat_ex, rho, theta)
% Outputs the first and second derivatives of the empirical log-likelihood 
% under a linear utility specification.
%
% g is the network, 
% stat_en is the vector of endogenous statistics (the output of en_stat.m), 
% stat_ex is the vector of exogenous statistics (the output of ex_stat.m), 
% theta_hat is the structural parameter.
%
% If there are no endogenous statistics, pass the empty matrix [].
%
% This function uses the multiprod function.

G = g(2:end,2:end) ;
N = size(g,1)-1; % number of agents
r0 = rho(1);    %Misspecified probabilities
r1 = rho(2);
p = size(stat_en,3) + size(stat_ex,3); % number of parameters 

c = - cat(2, eye(3), zeros(3,1))*([1- r0-r1 , 0 , 0, 0 ; 0, 1-r0-r1, 0, 0 ; 0 , 0, (1-r0-r1)^2, r0*(1-r0-r1); 0 , 0, 0 , 1-r0-r1 ]\[r0;r0;r0^2;r0] );
C =  cat(2, eye(3), zeros(3,1))/([1- r0-r1 , 0 , 0, 0 ; 0, 1-r0-r1, 0, 0 ; 0 , 0, (1-r0-r1)^2, r0*(1-r0-r1); 0 , 0, 0 , 1-r0-r1 ]);

% Computing Endogenous Statistics.
stat_en_sum = stat_en(:,:,3)+ stat_en(:,:,3)';
stat_en = cat(3, stat_en(:,:,1), stat_en(:,:,3), stat_en(:,:,2), stat_en_sum);
stat_en =  repmat( reshape(c,[1,1,3]), [N, N, 1] ) +  multiprod(stat_en, C', [2,3], [1,2]) ;

ES_ij = cat(3, stat_en, stat_ex);               % [N x N x p] ijth entry = \hat E[S_{ij} | X, sigma]
V_ij = multiprod(ES_ij, theta, [2 3], [1 2]);   % [N x N] ijth entry = \hat E[S_{ij} | X, \sigma]' \theta

% Likelihood
Phi     = normcdf(V_ij); 
Phi_mis =  r0 + (1-r0-r1).*Phi ; 
Phi_mis(Phi_mis==1)=0.999 ;
Phi_mis(Phi_mis==0)=0.001 ;
phi     = normpdf(V_ij);
phi(1:(size(phi,1)+1):end) = 0; % zero out diagonal elements

% Likelihood
L  =  - sum(sum(G.*log(Phi_mis) + (1-G).*log(1-Phi_mis),2),1)./(N*(N-1));

% First derivative of the log-likelihood.
dQ_summands = - multiprod(ES_ij, ( G - Phi_mis ).*(1-r0-r1).*phi./(Phi_mis.*(1-Phi_mis)) , 3) ; % [N x N x p]
dQ = squeeze(sum(sum(dQ_summands,1),2))/(N*(N-1)); % [p x 1]

