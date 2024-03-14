function T200_FullEst(rho, TE)

% Misclassified Probs.
Nr = size(rho, 2) ;
Theta_hat   = NaN( 19 , Nr);
SE          = NaN( 19 , Nr);

varnames = {'reciprocation'; ' '; 'in degree'; ' '; 'supported trust'; ' '; ...
            'constant'; ' '; 'same religion'; ' '; 'same sex'; ' '; ...
            'same caste'; ' '; 'same language'; ' '; 'same family'; ' '};


if TE==0
    for i=1:Nr
        load(['./results_mis/theta_hat_l0_r0_',num2str(rho(1,i)), '_r1_',num2str(rho(2,i)),'.mat'], 'theta_hat') 
        load(['./results_mis/Avar_hat_l0_r0_', num2str(rho(1,i)), '_r1_',num2str(rho(2,i)),'.mat'], 'V_hat')
        Theta(: , i)   = theta_hat         ;
        SE(: ,i)       = sqrt(diag(V_hat)) ;
    end
    
    headers_T =  {'$r_0=0$' ; '$r_0=0.01$'; '$r_0=0.02$'; '$r_0=0.03$' };
            
   
elseif TE==1
    
    for i=1:Nr
        load(['./results_mis/theta_hat_l0_r1_',num2str(rho(2,i)),'.mat'], 'theta_hat') 
        load(['./results_mis/Avar_hat_l0_r1_', num2str(rho(2,i)),'.mat'], 'V_hat')
        Theta(: , i)   = theta_hat         ;
        SE(: ,i)       = sqrt(diag(V_hat)) ;
    end
    
    headers_T =  {'$r_1=0$' ; '$r_1=0.05$'; '$r_1=0.1$' ; '$r_1=0.2$'; '$r_1=0.3$' ; ...
                '$r_1=0.4$'; '$r_1=0.5$'; '$r_1=0.6$'; '$r_1=0.7$' ; '$r_1=0.8$'; '$r_1=0.9$'};
    
else
    for i=1:Nr
        load(['./results_mis/theta_hat_l0_r0_',num2str(rho(1,i)), '_r1_',num2str(rho(2,i)),'.mat'], 'theta_hat') 
        load(['./results_mis/Avar_hat_l0_r0_', num2str(rho(1,i)), '_r1_',num2str(rho(2,i)),'.mat'], 'V_hat')
        Theta(: , i)   = theta_hat         ;
        SE(: ,i)       = sqrt(diag(V_hat)) ;
    end
    
    headers =  {'$r_0=0.01$, $r_1=0.05$'; '$r_0=0.01$, $r_1=0.1$'; '$r_0=0.01$, $r_1=0.2$$' ; ...
                '$r_0=0.02$, $r_1=0.05$' ; '$r_0=0.02$, $r_1=0.1$' ; '$r_0=0.02$, $r_1=0.2$'   ; ...
                '$r_0=0.03$, $r_1=0.05$' ; '$r_0=0.03$, $r_1=0.1$' ; '$r_0=0.03$, $r_1=0.2$' };
end

% Aggregate Statistics.
    End_gamma   = Theta(1:3, :) ;
    End_SE      = SE(1:3, :)    ;
    Hom_beta    = [Theta(4 , :); Theta(15:end , :)] ;
    Hom_SE      = [SE(4, :);  SE(15:end, :) ]       ;
    
    theta    =  cat(1, End_gamma,  Hom_beta) ;
    theta_SE =  cat(1, End_SE,  Hom_SE)      ;
    
    % Output Tables 
    output = NaN(size(theta,1)*2, Nr) ;     % Coefficients and SE.
    
for i= 1:size(theta,1) 

    output(2*i-1, :)   = theta(i,:)    ;
    output(2*i  , :)   = theta_SE(i,:) ;
    
end


if TE==0
    
   results =     output;
   headers =  headers_T;
    
    % Table1 : MLE estimates
    matrix2latex(results,'./results_mis/Table_MLE_r0.tex','rowLabels', varnames, 'columnLabels', headers, 'alignment', 'c','format', '%7.3f')
    
    
elseif TE==1
    
    results      = output(:,1:7);
    results_Ext  = cat(2, output(:,1), output(:,8:end) ) ;
    
    headers      = headers_T(1:7);     
    headers_Ext  = headers_T; headers_Ext(2:7)=[];
    
    % Table1 : MLE estimates
    matrix2latex(results,'./results_mis/Table_MLE_r1.tex','rowLabels', varnames, 'columnLabels', headers, 'alignment', 'c','format', '%7.3f')
    % Table2 : Extension
    matrix2latex(results_Ext,'./results_mis/Table_ExtMLE_r1.tex','rowLabels', varnames, 'columnLabels', headers_Ext, 'alignment', 'c','format', '%7.3f')
    
else

    matrix2latex(output,'./results_mis/Table_MLE_r0r1.tex','rowLabels', varnames, 'columnLabels', headers, 'alignment', 'c','format', '%7.3f')
    

end
    
%--------------------------------------------------------------------------
    % Print Table
    if TE==0||TE==1
    
    fprintf('Tables Generated, Printing Results\n')
    fprintf('Type Error =%.2g \n', TE)
        
    T_MLE      = array2table(results,     'VariableNames', headers) ;
    fprintf('Table: MLE Estimates\n')
    disp(T_MLE)
   
        
    else 
        fprintf('Table: Mixed Completed\n')
    end
%--------------------------------------------------------------------------


end

