function T300_ConfInt(rho, alpha, TE)

% Misclassified Probs.
Nr = size(rho, 2) ;
Theta   = NaN( 19 , Nr);
SE      = NaN( 19 , Nr);

varnames = {'reciprocation'; 'in degree'; 'supported trust'; ...
            'constant'; 'same religion'; 'same sex'; ...
            'same caste'; 'same language'; 'same family'};


if TE==0

    for i=1:Nr
        load(['./results_mis/theta_hat_l0_r0_',num2str(rho(1,i)), '_r1_',num2str(rho(2,i)),'.mat'], 'theta_hat') 
        load(['./results_mis/Avar_hat_l0_r0_', num2str(rho(1,i)), '_r1_',num2str(rho(2,i)),'.mat'], 'V_hat')
        Theta(: , i)   = theta_hat         ;
        SE(: ,i)       = sqrt(diag(V_hat)) ;
    end
    
    headers_T =  {'$LB(r_0=0)$'    ; '$UB(r_0=0)$'    ; '$LB(r_0<=0.01)$' ; '$UB(r_0<=0.01)$'; ...
                  '$LB(r_0<=0.02)$' ; '$UB(r_0<=0.02)$' ; '$LB(r_0<=0.03)$'  ; '$UB(r_0<=0.03)$' };
                
    headers_ratio_T  = {'$C(r_0=0.01)/C(r_0=0)$' ; '$C(r_0=0.02)/C(r_0=0)$' ; ...
                        '$C(r_0=0.03)/C(r_0=0)$' };

    
elseif TE==1
    
    for i=1:Nr
        load(['./results_mis/theta_hat_l0_r1_',num2str(rho(2,i)),'.mat'], 'theta_hat') 
        load(['./results_mis/Avar_hat_l0_r1_', num2str(rho(2,i)),'.mat'], 'V_hat')
        Theta(: , i)   = theta_hat         ;
        SE(: ,i)       = sqrt(diag(V_hat)) ;
    end
    
    headers_T ={'$LB(r_1=0)$'    ; '$UB(r_1=0)$'    ; '$LB(r_1<=0.05)$' ; '$UB(r_1<=0.05)$' ; ...
                '$LB(r_1<=0.1)$' ; '$UB(r_1<=0.1)$' ; '$LB(r_1<=0.2)$'  ; '$UB(r_1<=0.2)$'  ; ...
                '$LB(r_1<=0.3)$' ; '$UB(r_1<=0.3)$' ; '$LB(r_1<=0.4)$'  ; '$UB(r_0<=0.4)$'  ; ...
                '$LB(r_1<=0.5)$' ; '$UB(r_1<=0.5)$' ; '$LB(r_1<=0.6)$'  ; '$UB(r_1<=0.6)$'  ; ...
                '$LB(r_1<=0.7)$' ; '$UB(r_1<=0.7)$' ; '$LB(r_1<=0.8)$'  ; '$UB(r_1<=0.8)$'  ; ...
                '$LB(r_1<=0.9)$' ; '$UB(r_1<=0.9)$'};
            
    headers_ratio_T  = {'$C(r_1=0.05)/C(r_1=0)$' ; '$C(r_1=0.1)/C(r_1=0)$' ; ...
                        '$C(r_1=0.2)/C(r_1=0)$'  ; '$C(r_1=0.3)/C(r_1=0)$' ; ...
                        '$C(r_1=0.4)/C(r_1=0)$'  ; '$C(r_1=0.5)/C(r_1=0)$' ; ...
                        '$C(r_1=0.6)/C(r_1=0)$'  ; '$C(r_1=0.7)/C(r_1=0)$' ; ...
                        '$C(r_1=0.8)/C(r_1=0)$'  ; '$C(r_1=0.9)/C(r_1=0)$' };


else %mixed
        load(['./results_mis/theta_hat_l0_r0_',num2str(0), '_r1_',num2str(0),'.mat'], 'theta_hat') 
        load(['./results_mis/Avar_hat_l0_r0_', num2str(0), '_r1_',num2str(0),'.mat'], 'V_hat')
        Theta_base   = cat(1, theta_hat(1:4), theta_hat(15:end)  )    ;
        SE_base      = sqrt(diag(V_hat)) ;
        SE_base      = cat(1, SE_base(1:4), SE_base(15:end))  ;
        
    for i=1:Nr
        load(['./results_mis/theta_hat_l0_r0_',num2str(rho(1,i)), '_r1_',num2str(rho(2,i)),'.mat'], 'theta_hat') 
        load(['./results_mis/Avar_hat_l0_r0_', num2str(rho(1,i)), '_r1_',num2str(rho(2,i)),'.mat'], 'V_hat')
        Theta(: , i)   = theta_hat         ;
        SE(: ,i)       = sqrt(diag(V_hat)) ;
    end
    
headers=  {'$LB(r0=0.01, r1=0.05)$' ; '$UB(r0=0.01, r1=0.05)$'; '$LB(r0=0.01, r1=0.1)$' ; '$UB(r0=0.01, r1=0.1)$'; '$LB(r0=0.01, r1=0.2)$' ; '$UB(r0=0.01, r1=0.2)$'; '$LB(r0=0.01, r1=0.3)$' ; '$UB(r0=0.01, r1=0.3)$'; '$LB(r0=0.01, r1=0.4)$' ; '$UB(r0=0.01, r1=0.4)$'; ...
           '$LB(r0=0.02, r1=0.05)$' ; '$UB(r0=0.02, r1=0.05)$'; '$LB(r0=0.02, r1=0.1)$' ; '$UB(r0=0.02, r1=0.1)$'; '$LB(r0=0.02, r1=0.2)$' ; '$UB(r0=0.02, r1=0.2)$'; '$LB(r0=0.02, r1=0.3)$' ; '$UB(r0=0.02, r1=0.3)$'; '$LB(r0=0.02, r1=0.4)$' ; '$UB(r0=0.02, r1=0.4)$'; ...
           '$LB(r0=0.03, r1=0.05)$' ; '$UB(r0=0.03, r1=0.05)$'; '$LB(r0=0.03, r1=0.1)$' ; '$UB(r0=0.03, r1=0.1)$'; '$LB(r0=0.03, r1=0.2)$' ; '$UB(r0=0.03, r1=0.2)$'; '$LB(r0=0.03, r1=0.3)$' ; '$UB(r0=0.03, r1=0.3)$'; '$LB(r0=0.03, r1=0.4)$' ; '$UB(r0=0.03, r1=0.4)$'};
                
    
headers_ratio  = {'$C(r0=0.01, r1=0.05)/C(r0=0.0, r1=0.0)$'  ; '$C(r0=0.01, r_1=0.1)/C(r0=0.0, r1=0.0)$' ; ...
                    '$C(r0=0.01, r_1=0.2)/C(r0=0.0, r1=0.0)$'  ; '$C(r0=0.01, r_1=0.3)/C(r0=0.0, r1=0.0)$' ; ...
                    '$C(r0=0.01, r_1=0.4)/C(r0=0.0, r1=0.0)$'  ; ...
                    '$C(r0=0.02, r1=0.05)/C(r0=0.0, r1=0.0)$'  ; '$C(r0=0.02, r_1=0.1)/C(r0=0.0, r1=0.0)$' ; ...
                    '$C(r0=0.02, r_1=0.2)/C(r0=0.0, r1=0.0)$'  ; '$C(r0=0.02, r_1=0.3)/C(r0=0.0, r1=0.0)$' ; ...
                    '$C(r0=0.02, r_1=0.4)/C(r0=0.0, r1=0.0)$'  ; ...
                    '$C(r0=0.03, r1=0.05)/C(r0=0.0, r1=0.0)$'  ; '$C(r0=0.03, r_1=0.1)/C(r0=0.0, r1=0.0)$' ; ...
                    '$C(r0=0.03, r_1=0.2)/C(r0=0.0, r1=0.0)$'  ; '$C(r0=0.03, r_1=0.3)/C(r0=0.0, r1=0.0)$' ; ...
                    '$C(r0=0.03, r_1=0.4)/C(r0=0.0, r1=0.0)$'  };
                
end

End_gamma   = Theta(1:3, :) ;
End_SE      = SE(1:3, :)    ;
Hom_beta    = [Theta(4 , :); Theta(15:end , :)] ;
Hom_SE      = [SE(4, :);  SE(15:end, :) ]       ;

theta    =  cat(1, End_gamma,  Hom_beta) ;
theta_SE =  cat(1, End_SE,  Hom_SE)      ;

if TE == 2
    theta    = cat(2, Theta_base, theta) ; 
    theta_SE = cat(2, SE_base, theta_SE)  ;
else 
fprintf('Continue \n')
end

% 1/ Confidence Intervals
quant = norminv(1-alpha/2, 0, 1) ;
LB = theta - theta_SE.*quant ; 
UB = theta + theta_SE.*quant ; 

% 2/ Union of Confidence Intervals. 
LB_Union = NaN(size(theta,1), size(theta,2))     ; 
UB_Union = NaN(size(theta,1), size(theta,2))     ; 
CI_Union = NaN(size(theta,1), 2*size(theta,2))   ; 

LB_Union(:,1) =  LB(:,1);
UB_Union(:,1) =  UB(:,1);
CI_Union(:,1) =  LB_Union(:,1);
CI_Union(:,2) =  UB_Union(:,1);

for i=2:size(theta,2)
    LB_Union(:,i)       = min(LB(:,1:i),[], 2);
    UB_Union(:,i)       = max(UB(:,1:i),[], 2);
    CI_Union(:,2*i-1)   = LB_Union(:,i) ;
    CI_Union(:,2*i)     = UB_Union(:,i) ;
end

% 3/ Ratios
ratio_lengths = (UB_Union-LB_Union)./(UB(:,1)-LB(:,1)); 
ratio_lengths(:,1) =[];

    if TE==0

    % Store results for tables
    results         = CI_Union;
    results_ratios  = ratio_lengths;
            
    headers      = headers_T;     
    headers_ratio = headers_ratio_T;
 
    
    % Table 1. Union of Confidence Intervals
    matrix2latex(results, './results_mis/Table_CIsUnion_r0.tex','rowLabels', varnames, 'columnLabels', headers, 'alignment', 'c','format', '%7.3f')
    % Table 2. Ratio of Lengths
    matrix2latex(results_ratios, './results_mis/Table_CIsRatios_r0.tex','rowLabels', varnames, 'columnLabels', headers_ratio, 'alignment', 'c','format', '%7.3f')

    elseif TE==1
        
    % Store results for tables
    results         = CI_Union(:,1:14);
    results_ratios  = ratio_lengths(:,1:6);

    results_Ext         = cat(2, CI_Union(:,1:2), CI_Union(:,15:end) ) ;
    results_ratios_Ext  = ratio_lengths(:,7:end);
        
    headers      = headers_T(1:14);     
    headers_Ext  = headers_T; headers_Ext(3:14)=[];    

    headers_ratio = headers_ratio_T(1:6);
    headers_ratio_Ext = headers_ratio_T;  headers_ratio_Ext(1:6)=[];
    
    % Table 1. Union of Confidence Intervals
    matrix2latex(results, './results_mis/Table_CIsUnion_r1.tex','rowLabels', varnames, 'columnLabels', headers, 'alignment', 'c','format', '%7.3f')
    % Table 2. Ratio of Lengths
    matrix2latex(results_ratios, './results_mis/Table_CIsRatios_r1.tex','rowLabels', varnames, 'columnLabels', headers_ratio, 'alignment', 'c','format', '%7.3f')
    % Table 3: Extension Unions
    matrix2latex(results_Ext, './results_mis/Table_ExtCIsUnion_r1.tex','rowLabels', varnames, 'columnLabels', headers_Ext, 'alignment', 'c','format', '%7.3f')
    % Table 4. Extension Ratios
    matrix2latex(results_ratios_Ext, './results_mis/Table_ExtCIsRatios_r1.tex','rowLabels', varnames, 'columnLabels', headers_ratio_Ext, 'alignment', 'c','format', '%7.3f')
    
    else % Mixed
    
    CI_Union(:,1:2) = [];
    results_01 = CI_Union(:, 1:10)   ; 
    results_02 = CI_Union(:, 11:20)  ; 
    results_03 = CI_Union(:, 21:30)  ; 

    results_ratios_01   = ratio_lengths(:,1:5)  ;  
    results_ratios_02   = ratio_lengths(:,6:10) ;  
    results_ratios_03   = ratio_lengths(:,11:15);  
                
    % Table 1. Union of Confidence Intervals
    matrix2latex(results_01, './results_mis/Table_CIsUnion_r0r1_001.tex','rowLabels', varnames, 'columnLabels', headers(1:10), 'alignment', 'c','format', '%7.3f')
    matrix2latex(results_02, './results_mis/Table_CIsUnion_r0r1_002.tex','rowLabels', varnames, 'columnLabels', headers(11:20), 'alignment', 'c','format', '%7.3f')
    matrix2latex(results_03, './results_mis/Table_CIsUnion_r0r1_003.tex','rowLabels', varnames, 'columnLabels', headers(21:30), 'alignment', 'c','format', '%7.3f')
    % Table 2. Ratio of Leng
    matrix2latex(results_ratios_01, './results_mis/Table_CIsRatios_r0r1_001.tex','rowLabels', varnames, 'columnLabels', headers_ratio(1:5), 'alignment', 'c','format', '%7.3f')
    matrix2latex(results_ratios_02, './results_mis/Table_CIsRatios_r0r1_002.tex','rowLabels', varnames, 'columnLabels', headers_ratio(6:10), 'alignment', 'c','format', '%7.3f')
    matrix2latex(results_ratios_03, './results_mis/Table_CIsRatios_r0r1_003.tex','rowLabels', varnames, 'columnLabels', headers_ratio(11:15), 'alignment', 'c','format', '%7.3f')
    
    end
    
    % Print Table
    if TE==0||TE==1
    
    fprintf('Tables Generated, Printing Results\n')
    fprintf('Type Error =%.2g \n', TE)
        
    T_Union      = array2table(results, 'VariableNames', headers) ;
    T_Ratios     = array2table(results_ratios, 'VariableNames', headers_ratio);
    fprintf('Table: Union Confidence Intervals\n')
    disp(T_Union) 
    fprintf('Table: Ratios of Lenghts \n')
    disp(T_Ratios)
   
   
        
    else 
        fprintf('Table: Mixed Completed\n')
    end
    
end

