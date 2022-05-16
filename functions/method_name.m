function n = method_name(m)

% old names
% -----------------------------------
    if strcmp(m,'res_no_reg')
        
        n = 'cob';
        
    elseif strcmp(m,'no_reg_rebal')
        
        n = 'cob + rebal';
        
    elseif strcmp(m,'reg_only_rebal')
        
        n = 'bal + tr inv + rebal';        
        
%     elseif strcmp(m,'tr_inv_rebal')
%         
%         n = 'cob + tr inv + rebal';          
        
    elseif strcmp(m,'tr_inv_only_rebal')
        
        n = 'tr inv + rebal';        
        
% new names
% -----------------------------------   
    elseif strcmp(m,'cob')
        
        n = 'cob';  
        
    elseif strcmp(m,'cob_rebal')
        
        n = 'cob + rebal';   
        
    elseif strcmp(m,'tr_inv_rebal')
        
        n = 'tr inv + rebal';  
        
    elseif strcmp(m,'cob_tr_inv_rebal')
        
        n = 'cob + tr inv + rebal';    
        
    elseif strcmp(m,'cob_tr_inv_rebal_new_lambda')
        
        n = 'cob + tr inv + rebal (new lambda)';  
        
    elseif strcmp(m,'cob_tr_inv_rebal_var_lambda')
        
        n = 'cob + tr inv + rebal (variable lambda)';  
        
    elseif strcmp(m,'no_cob_rebal_analysis')
        
        n = 'cob + rebal';   
        
    elseif strcmp(m,'cob_rebal_analysis')
        
        n = 'cob + rebal';   
        
    elseif strcmp(m,'cob_rebal_analysis_extra')
        
        n = 'cob + rebal';           
        
    elseif strcmp(m,'cob_rebal_analysis_ls')
        
        n = 'cob + rebal (ls)';           
        
    elseif strcmp(m,'cob_tr_inv_rebal_var_lambda_analysis')
        
        n = 'cob + tr inv + rebal (variable lambda)';        
        
    elseif strcmp(m,'cob_tr_inv_rebal_var_lambda_analysis_extra')
        
        n = 'cob + tr inv + rebal (variable lambda)';  
        
% neurips names
% -----------------------------------  

    elseif strcmp(m,'gd')
        
        n = 'gd';  
        
    elseif strcmp(m,'gdr')
        
        n = 'gdr';  
        
    elseif strcmp(m,'irpg')
        
        n = 'irpg';          
        
    else
        n = 'unknown';
    end



end

