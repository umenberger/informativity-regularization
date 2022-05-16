function n = method_name_paper(m)

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
        
        n = '$\mathcal{L}_\mathtt{OE}$';  
        
    elseif strcmp(m,'cob_rebal')
        
        n = '$\mathcal{L}_\mathtt{OE}$ w/ recond.';   
        
    elseif strcmp(m,'tr_inv_rebal')
        
        n = 'tr inv + rebal';  
        
    elseif strcmp(m,'cob_tr_inv_rebal')
        
        n = 'cob + tr inv + rebal';    
        
    elseif strcmp(m,'cob_tr_inv_rebal_new_lambda')
        
        n = 'cob + tr inv + rebal (new lambda)';  
        
    elseif strcmp(m,'cob_tr_inv_rebal_var_lambda')
        
        n = '$\mathcal{L}_\mathtt{OE} + \lambda \mathcal{R}_\mathtt{info}$ w/ recond.';  
        
% neurips names
% -----------------------------------  

    elseif strcmp(m,'gd')
        
        n = '$\mathcal{L}_\mathtt{OE}$';  
        
    elseif strcmp(m,'gdr')
        
        n = '$\mathcal{L}_\mathtt{OE}$ w/ recond.';  
        
    elseif strcmp(m,'irpg')
        
        n = '$\mathcal{L}_\mathtt{OE} + \lambda \mathcal{R}_\mathtt{info}$ w/ recond.';         
        
    else
        n = 'unknown';
    end



end

