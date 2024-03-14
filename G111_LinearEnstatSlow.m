function f = G111_LinearEnstatSlow(g, Xest, Z, lambda)
% Same as function en_stat, except certain memory-intensive vectorized
% calculations are replaced with slow for-loops.

N = size(g,1) - 1;
fprintf('    en_stat_slow, N=%d\n',N);

tempX = dataset({Xest,'pid','drop','age','gender','HOH','hindu','caste2','caste3','educ_primary','educ_secondary','educ_puc','educ_ideg','educ_deg','educ_oth','hindi','kannada','malayalam','marati','tamil','telugu','urdu','english'});
tempID = dataset({g(2:(N+1))','pid'});
Xg = double(join(tempID, tempX, 'Type', 'leftouter'));
%Xg = Xg(:,3:size(Xg,2));% omit pid, drop
Xg = Xg(:,[4:8 15:size(Xg,2)]);% omit pid, drop, age, educ
Z = Z(2:size(Z,1),2:size(Z,2)); % strip pids

%% Create vectors of endogenous statistics
% [N^2, N] where N^2 captures the [Ekikj]i,j=<N for a given k=1,..,N:  
disp('        Create endogenous stats.');
G = g(2:(N+1),2:(N+1)); % strip away pids

stat_3 = reshape ( repmat(permute(G,[3 2 1]), [N 1 1]) .* repmat( permute(G,[2 3 1]), [1 N 1]), [N^2 N] );
stat_4 = reshape ( repmat(permute(G,[3 2 1]), [N 1 1]) + repmat( permute(G,[2 3 1]), [1 N 1]), [N^2 N] );


%% Create first-step estimators
% Attributes
categ = cat(3, repmat(permute(Xg(:,2:size(Xg,2)),[1 3 2]),[1 N 1]), repmat(permute(Xg(:,2:size(Xg,2)),[3 1 2]),[N 1 1])); % [N x N x (d-1)*2] where d is the dimension of X_1. ijth element is the categorical components of vector (X_i,X_j)
% ID.
ordrd = cat(3, repmat(permute(Xg(:,1),[1 3 2]),[1 N 1]), repmat(permute(Xg(:,1),[3 1 2]),[N 1 1])); % [N x N x 2] ijth element is (age_i, age_j), the only ordered components of (X_i,X_j)
clear tempX tempID Xg 

recip_lnr     = zeros(N,N,N);
indeg_lnr     = zeros(N,N,N);
supptrust_lnr = zeros(N,N,N);
doubindeg_lnr = zeros(N,N,N);

parfor k=1:N %Main looping along 3rd dim.

for i=1:N
    for j=1:N
        if(i==j)
            continue;
        end
        
        weight = lambda .^ (sum(categ ~= repmat(categ(i,j,:),[N N 1]),3) + sum(abs(ordrd - repmat(ordrd(i,j,:),[N N 1])),3) + sum(Z ~= repmat(Z(i,j,:),[N N 1]),3)); % [N x N] ijth entry = number of differences between categorical components of (X_k,X_l) and (X_i,X_j) + |age_i - age_k| + |age_j - age_l| + 1{Z_ij ~= Z_kl}
        denom  = sum(sum(weight,1),2)/N^2 ;
        
        recip_lnr(i,j,k)         =   (1/N)*G(k,:)*weight(:,k)/denom   ;
        indeg_lnr(i,j,k)         =   (1/N^2)*sum(weight* G(k, :)')/denom   ; %sum(weight* G(k,:)')
        supptrust_lnr(i,j,k)     =   (1/N^2)* sum(sum(weight.* reshape(stat_3(:,k), [N N ]) ))/denom  ;
        doubindeg_lnr(i,j,k)     =   (1/N^2)* sum(sum(weight.* reshape(stat_4(:,k), [N N ]) ))/denom  ;
        
    end
end

end 

recip_lnr       =  reshape(recip_lnr, [N^2 N] ) ;
indeg_lnr       =  reshape(indeg_lnr, [N^2 N] ) ;
supptrust_lnr   =  reshape(supptrust_lnr, [N^2 N] ) ;
doubindeg_lnr   =  reshape(doubindeg_lnr, [N^2 N] ) ;

f = cat(3, recip_lnr, indeg_lnr, supptrust_lnr, doubindeg_lnr);

