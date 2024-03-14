function T100_Descriptive(vils)

% Summary statistics

load('X.mat');

%vils = [6 12 29 34 35 46 71 74 76];

stats      = zeros(length(vils),6);          %Covariates 
in_degree  = zeros(303, length(vils));      
out_degree = zeros(303, length(vils));   
Nv = zeros(length(vils),1) ;

for i=1:length(vils)
    
    g = csvread(['directed_adjacency_matrices/lendmoney',num2str(vils(i)),'.csv']);
    G = g(2:end, 2:end) ;
    N = size(G,1); Nv(i) = N ;
    
    stats(i,1) = N;  %Number of individuals.
    
    Xs = zeros( N, 5);
    for k=1:N
        Xs(k,:) = X(X(:,1)==g(k+1,1),[2:3 5 9:10]);
    end
    stats(i,2) = mean(Xs(:,1));     % average age
    stats(i,3) = mean(Xs(:,2));     % pct female
    stats(i,4) = mean(Xs(:,3));     % pct hindu
    stats(i,5) = mean(Xs(:,4));     % pct OBC
    stats(i,6) = mean(Xs(:,5));     % pct sched
        
    ideg  = sum(G,1)' ; %In-degree 
    odeg  = sum(G,2)  ; %Out-degree 
    
    slack = zeros(303-N, 1) ;
    
    in_degree(:, i)  = cat(1, ideg, slack) ;
    out_degree(:, i) = cat(1, odeg, slack) ;    
   
end

Ideg = in_degree(:); Odeg = out_degree(:);
mean_id  = sum(Ideg)/sum(Nv) ; 
var_id   = sum((Ideg - mean_id).^2)/(2031-1);
mean_od  = sum(Odeg)/sum(Nv) ; 
var_od   = sum((Odeg - mean_od).^2)/(2031-1);

stat_in  = [ mean_id, sqrt(var_id),  min(Ideg), max(Ideg)] ;
stat_out = [ mean_od, sqrt(var_od),  min(Odeg), max(Odeg)] ;

res_stats  = [mean(stats)', std(stats)', min(stats)', max(stats)'];

% Export tables to LaTex.
results = cat(1, res_stats(1,:), [stat_in; stat_out], res_stats(2:end,:) ) ;

headers =  {'mean' ; 'sd' ; 'min'; 'max' };
    
varnames = {'villagers'; 'in_degree'; 'out_degree'; 'average age'; 'share female';  'share Hindu';  'share OBC';  'share scheduled'};

matrix2latex(results,'./results_mis/Table_descriptive.tex','rowLabels', varnames, 'columnLabels', headers, 'alignment', 'c','format', '%7.3f')

%--------------------------------------------------------------------------
% Print Table

fprintf('Tables Generated, Printing Results\n')

T_stats  = array2table(results,     'VariableNames', headers) ;

fprintf('Table: Descriptive Statistics\n')
disp(T_stats)



end

