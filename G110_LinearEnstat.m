function f = G110_LinearEnstat(g, Xest, Z, lambda)
% Outputs a [N^2 x N x 4] matrix. 
% Third dimension: linear representation of 4 network statistics.
% Second dimension: Linear representation for each k=1, ..., N.
% First Dimension: (i,j)-specific value.

% Xest is the matrix created in format_covariates.
% Z is the family adjacency matrix, including person ids (pids).
% g is the directed adjacency matrix.
% lambda is the smoothing parameter.
%
% Works for villages under 220 people with around 10gb memory.

%load('X.mat');
%Z = csvread(['directed_adjacency_matrices/rel',num2str(w),'.csv']);

N = size(g,1) - 1;
fprintf('    en_stat, N=%d\n',N);

tempX   = dataset({Xest,'pid','drop','age','gender','HOH','hindu','caste2','caste3','educ_primary','educ_secondary','educ_puc','educ_ideg','educ_deg','educ_oth','hindi','kannada','malayalam','marati','tamil','telugu','urdu','english'});
tempID  = dataset({g(2:(N+1))','pid'});
Xg      = double(join(tempID, tempX, 'Type', 'leftouter'));
Xg      = Xg(:,[4:8 15:size(Xg,2)]); % omit pid, drop, age, educ: [Nx14]
Z       = Z(2:size(Z,1),2:size(Z,2)); % strip pids: [NxN]

%% Create weights

icat     = repmat( permute( Xg(:,2:size(Xg,2)),[1 3 2]) , [1 N 1]); % categorical attributes of i
jcat     = repmat( permute( Xg(:,2:size(Xg,2)), [3 1 2]), [N 1 1]); % categorical attributes of j
cat_diff = sum(icat~=jcat,3); % [N x N] number of disagreeing categorical components between X_i, X_j
iord     = repmat(permute(Xg(:,1),[1 3 2]),[1 N 1]); % ordered attributes of i
jord     = repmat(permute(Xg(:,1),[3 1 2]),[N 1 1]); % ordered attributes of j
ord_diff = sum(abs(iord-jord),3); % [N x N] sum of absolute differences between ordered components of X_i, X_j

fprintf('        Create diffs_X. ');
tstart=tic;
q = 1:N; 
q = q(ones(1,N),:);
tmp = cat_diff + ord_diff;
clear tempX tempID icat jcat cat_diff iord jord ord_diff 

Xdiff_tmp   = repmat( tmp(:) , [1 N^2] );
diffs_X     = Xdiff_tmp + Xdiff_tmp';
    % [N^2 x N^2] first memory bottleneck.
    % (ij,kl)th entry:
    % First component of the sum is number of disagreeing components between the 
    % categorical components of X_i and X_k, plus the absolute difference
    % between the ordered components.
    % Second component is the same thing for X_j and X_l.
clear tmp q Xdiff_tmp
telapsed = toc(tstart);
fprintf('Time elapsed: %.3g.\n', telapsed);

fprintf('        Create diffs. ');
tstart=tic;
tmp = repmat(Z(:),[1 size(Z,1)^2]);
diffs_fam = tmp ~= tmp'; % [N^2 x N^2] matrix of indicators for whether or not Z_{ij}, Z_{kl} disagree

diffs = diffs_X + diffs_fam;
clear diffs_X diffs_fam tmp
telapsed = toc(tstart);
fprintf('Time elapsed: %.3g.\n', telapsed);


fprintf('        Create weight. ');
tstart=tic;
weight = lambda .^ diffs;  % w = lambda^d(Xik, Xil), [N^2, N^2]
%clear diffs
telapsed = toc(tstart);
fprintf('Time elapsed: %.3g.\n', telapsed);


%% Create Vector of Influence Functions

G = g(2:(N+1),2:(N+1));

stat_3 = reshape ( repmat(permute(G,[3 2 1]), [N 1 1]) .* repmat( permute(G,[2 3 1]), [1 N 1]), [N^2 N] );
        % [N^2 N]
stat_4 = reshape ( repmat(permute(G,[3 2 1]), [N 1 1]) +  repmat( permute(G,[2 3 1]), [1 N 1]), [N^2 N] );
        % [N^2 N]
          
%% Create weigted endogenous statistics.

disp('        Create endogenous stats.');

denom = (ones(1,N^2) * weight)'./N^2 ;    %[N^2 x 1]
denom(denom==0)=0.001;

recip_lnr     = zeros(N^2,N);
indeg_lnr     = zeros(N^2,N);
supptrust_lnr = zeros(N^2,N);
doubindeg_lnr = zeros(N^2,N);

for k=1:N
      recip_lnr(:,k)     = (( weight(:,((k-1)*N +1):k*N)*G(k,:)' )./denom)/N ;
      indeg_lnr(:,k)     = (( weight*reshape( repmat( G(k,:), [N,1]), [N^2, 1]) )./denom)/N^2;
      supptrust_lnr(:,k) = (( weight*stat_3(:,k) )./denom)/N^2  ;
      doubindeg_lnr(:,k) = (( weight*stat_4(:,k) )./denom)/N^2  ; 
end


f = cat(3, recip_lnr, indeg_lnr, supptrust_lnr, doubindeg_lnr);


