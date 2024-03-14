% Main m file for "Identification and Inference of Network Formation Games
% with Misclassified Links" by Luis E. Candelaria and Takuya Ura.
% Version: 10 March, 2022.
% 1. Initial setup 
%--------------------------------------------------------------------------
clear all
clc
format long

cd('./') %Set PATH to actual directory of this folder. 
addpath('code/Multiprod_2009');
addpath('code/matrix2latex');

vils = [6 12 29 34 35 46 71 74 76];
lambda = 0; L = 0;
TE=2;
% Misclassification probabilities

if TE==0
    
   r0 = [0, 0.01, 0.02, 0.03] ; 
   r1 = 0;
   rho = cat(1, r0, repmat(r1, [1, size(r0,2)]));
    
elseif TE==1
    
   r0 = 0 ;
   r1 = [0, 0.05, 0.1,  0.2,  0.3,  0.4,  0.5, 0.6, 0.7, 0.8, 0.9] ; 
   rho = cat(1, repmat(r0, [1, size(r1,2)]), r1);

else

  r0 = [0.01, 0.02,  0.03] ;
  r1 = [0.05, 0.1,  0.2,  0.3,  0.4];

  R0 = repmat(r0 , [size(r1,2),1]);
  R1 = repmat(r1' , [1,size(r0,2)]);

  rho = cat(2, R0(:), R1(:))' ; 


end
 
fprintf('Type Error=%.2g \n', TE)
fprintf('Distribution misclassification\n') 
disp(rho)

% Number of Parameters (3 endogenous and 16 exogenous)
load(['dMU/exstat_vil',num2str(vils(1)),'.mat']);
load(['dMU/enstat_vil',num2str(vils(1)),'_lam',num2str(lambda(1)),'.mat']);
p = size(stat_en,3) + size(stat_ex,3) ; 


%% 2. Create Network and Covariates
% Data Manipulation (Index D100_.m)
%--------------------------------------------------------------------------
 D100_gendirected(vils)      % generate adjacency matrices
 D100_format_covariates();   % generate matrix of individual characteristics
 D100_drop_notype_all(vils); % drop individuals with no characteristics (weren't surveyed) or missing crucial characteristics
% 
% Network and Homophily Covarites (Index D200_.m)
 D200_dMU([6 12 34 35 74], lambda, 0);   % villages with < 220 people
 D200_dMU([29 46 71 76], lambda, 1);     

%% 3. Create Asymptotic Linear Representation First Stage

G100_PsiGamma([6 12 34 35 74],lambda, 0) ;
G100_PsiGamma([29 46 71 76],lambda,  1) ;   % if insufficient memory, change '0' to '1' to use en_stat_slow instead


%% 4. Estimation of Parameters: 
    % E200_Qvillmis - Villages Aggregator.
    % E210_dQmis    - Computes Loss function and Gradient
%--------------------------------------------------------------------------
disp('Estimation parameters.');
options = optimoptions('fminunc', 'Display','iter', 'Algorithm','trust-region','SpecifyObjectiveGradient',true);

theta0 = zeros(p,1);

for r=1:size(rho,2) 
    
fprintf('r0=%.2g , r1=%.2g \n', rho(1, r), rho(2, r))

[theta_hat] = fminunc( @(theta) E200_Qvillmis(@E210_dQmis, vils, rho(:,r), theta), theta0, options) ;
    
save(['results_mis/theta_hat_l0_r0_',num2str(rho(1,r)),'_r1_',num2str(rho(2,r)) ,'.mat'], 'theta_hat');

end

%% 5. Asymptotic variance
% V100_Avarmis % V110_Avarnetmis

disp('Computing asymptotic variances.');
for r=1:size(rho,2) 
    fprintf('r0=%.2g , r1=%.2g \n', rho(1, r), rho(2, r))
    [V_hat]=V100_Avarmis(@V110_Avarnetmis,vils, rho(:,r));
end


%% 6. Tables Paper: False Negative
alpha =5/100;

% Descriptive Stats.
T100_Descriptive(vils)         

%Structural Parameters 
T200_FullEst(rho, TE)    

% Table Union of Confidence Intervals and Ratios.
T300_ConfInt(rho, alpha, TE )   

           

