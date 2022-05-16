function res = gd_filtering_state_cost(cntrl_init, sys, ops)

% This function performs gradient descent on the output error cost function
% Apologies: the code is not especially well-commented.
% Please write to umnbrgr@mit.edu for any clarification.

% Inputs:
%   cntrl_init - initial filter from which to begin local search
%   sys - true systems
%   ops - options for the search 

% Outputs:
%   res - a large matlab struct containing many results


% Extract search options from "ops" input
    if isfield(ops,'obj_tol')
        obj_tol = ops.obj_tol;
    else
        obj_tol = 1e-13;
        ops.obj_tol = obj_tol;
    end
    
    if isfield(ops,'grad_tol')
        grad_tol = ops.grad_tol;
    else
        grad_tol = 1e-6;
        ops.grad_tol = grad_tol;
    end    
    
    if isfield(ops,'termination_criterion')
        termination_criterion = ops.termination_criterion;
    else
        termination_criterion = 'gradient';
        ops.termination_criterion = termination_criterion;
    end    
    
    if strcmp(termination_criterion,'objective')
        term_criterion = 1;
    elseif strcmp(termination_criterion,'gradient')
        term_criterion = 2;
    else
        error('no such term criterion')
    end
                    
    if isfield(ops,'linesearch')
        linesearch = ops.linesearch;
    else
        linesearch = 1;
        ops.linesearch = linesearch;
    end   
    
    if isfield(ops,'num_iters')
        num_iters = ops.num_iters;
    else
        num_iters = 1e5;
        ops.num_iters = num_iters;
    end   
    
    if isfield(ops,'reg_type')
        reg_type = ops.reg_type;
    else
        reg_type = 0;
        ops.reg_type = reg_type;
    end   
    
    if isfield(ops,'lambda')
        lambda = ops.lambda;
    else
        lambda = 0.1;
        ops.lambda = lambda;
    end  
    
    lambda_c = 1;
    
    if isfield(ops,'norm_init_X22')
        norm_init_X22 = ops.norm_init_X22;
    else
        norm_init_X22 = 0;
        ops.norm_init_X22 = norm_init_X22;
    end     
    
    if isfield(ops,'rebalance_X22')
        rebalance_X22 = ops.rebalance_X22;
    else
        rebalance_X22 = 0;
        ops.rebalance_X22 = rebalance_X22;
    end      
    
    if isfield(ops,'save_progress')
        save_progress = ops.save_progress;
    else
        save_progress = 0;
        ops.save_progress = save_progress;
    end    
    
    if isfield(ops,'analysis')
        analysis = ops.analysis;
    else
        analysis = 0;
        ops.analysis = analysis;
    end      
    
    report_frequency_sec = 2*60;
    last_report_time = tic;

% Name of data file if "save progress" is used   
    dt = datestr(now,'yyyy_mm_dd_HH_MM_SS');
    fname = [dt '_inc_save'];    
     
    [ny,nx] = size(sys.C);
    nc = size(cntrl_init.A,1);
    
% change of basis to ensure X22 = I from first iteration.    
    if norm_init_X22
       
        [~,~,X] = kf_state_cost_ct(ss2mat(cntrl_init),sys);
        
        L = chol(X(nx+1:end,nx+1:end));
        L = L';

        cntrl_init_ = cntrl_init;
        cntrl_init_.A = L\cntrl_init.A*L;
        cntrl_init_.B = L\cntrl_init.B;
        cntrl_init_.C = cntrl_init.C*L;

        cntrl_init = cntrl_init_;
        
        fprintf('Using normalized X22\n')
                     
    end
    
    % initialization
    cntrl = ss2mat(cntrl_init);    

    cs = nan(1,num_iters); % output error cost
    fs = nan(1,num_iters); % regularized cost
    
    Rcs = nan(1,num_iters);
    dRc_norms = nan(1,num_iters);
    
    R12s = nan(1,num_iters);
    dR12_norms = nan(1,num_iters);  
    
    if reg_type == 5
        Ris = nan(1,num_iters);
        dRi_norms = nan(1,num_iters);
    end    
    
    min_sv_X12 = nan(1,num_iters);
    min_eig_X22 = nan(1,num_iters);
    error_eye_X22 = nan(1,num_iters);

    grad_c_norms = nan(1,num_iters); % gradient of cost
    grad_f_norms = nan(1,num_iters); % gradient of regularized cost
    
    grad_c_fro_norms = nan(1,num_iters);

    ctrbs = nan(1,num_iters);
    obsvs = nan(1,num_iters);
    
    stepsizes = nan(1,num_iters);
    
    if analysis
       
        [ny,nx] = size(sys.C);
        [cntrl_opt,~,~] = kf_via_care(sys.A,sys.C,sys.W,sys.V);
        cntrl_opt.C = eye(nx);
        cntrl_opt.D = zeros(nx,ny);   
        
        [c_opt, dc_opt, Sigma_opt] = kf_state_cost_ct(ss2mat(cntrl_opt),sys); 
        
        [~,~,~,~,~,nu_e_opt] = extended_Phi_map(cntrl_opt.A,cntrl_opt.B,cntrl_opt.C,Sigma_opt,sys);
                
        param_norms = nan(1,num_iters);
        param_norms_cvx = nan(1,num_iters);        
        param_error_cvx = nan(1,num_iters);
        min_eig_hessian = nan(1,num_iters);
        max_eig_hessian = nan(1,num_iters);
        min_svd_dPhi = nan(1,num_iters);
        max_svd_dPhi = nan(1,num_iters);  
        
        Kcvx_fro_norms = nan(1,num_iters);  
        Lcvx_fro_norms = nan(1,num_iters);  
        Mcvx_fro_norms = nan(1,num_iters);  
        Xcvx_fro_norms = nan(1,num_iters);  
        Ycvx_fro_norms = nan(1,num_iters); 
        
        Ucvx_2_norms = nan(1,num_iters); 
        Vcvx_2_norms = nan(1,num_iters); 
        
    end

    iter = 1;
    run = 1;
%     crossed_over = 0;
%     crossed_index = nan;
    term_msg = 'max iters exceeded';

    if linesearch

        LS_gamma = 0.2;
        LS_c = 0.01;
        
    else
        
        step_size = 0.1;

    end


    % multiple possible "regularized" objectives to optimize 

    if reg_type == 0 

        target = @(x) kf_state_cost_ct(x,sys);

    elseif reg_type == 1 

        target = @(x) kf_state_cost_ct(x,sys) + lambda*cntrl_state_reg(x,sys);  

    elseif reg_type == 2

        target = @(x) kf_state_cost_ct(x,sys) + lambda*tr_inv_ct_kz(x,sys);

    elseif reg_type == 3

        target = @(x) kf_state_cost_ct(x,sys) + lambda*(lambda_c*cntrl_state_reg(x,sys) + x_xinv_x_inv_reg(x,sys));

    elseif reg_type == 4

        target = @(x) cntrl_state_reg(x,sys) + x_xinv_x_inv_reg(x,sys);
        
    elseif reg_type == 5

        target = @(x) kf_state_cost_ct(x,sys) + lambda*(state_cov_eye(x,sys) + x_xinv_x_inv_reg(x,sys));
        
    elseif reg_type == 6
        
        target = @(x) tr_inv_ct_kz(x,sys);   
        
    elseif reg_type == 7
        
        target = @(x,l) kf_state_cost_ct(x,sys) + l*tr_inv_ct_kz(x,sys);
                
    else 
        error('no such reg')
    end
    
    lsp_counter = 0;

    
% start the local search
% -------------------------------------------------------------------------
    
    try
           
        while run
            
        % rebalance X22 if desired
            if rebalance_X22
               
                cntrl = ss2mat(balance_X22(mat2ss(cntrl,nx),sys));
                
            end            
            
            cntrl_sys = mat2ss(cntrl,nx);
            
        % compute controllability/observability measures for the current iterate
            ctrbs(iter) = min(svd(ctrb(cntrl_sys.A,cntrl_sys.B)));
            obsvs(iter) = min(svd(obsv(cntrl_sys.A,cntrl_sys.C)));

        % compute the cost 
            [c, dc, X] = kf_state_cost_ct(cntrl,sys);
            X12 = X(1:nx,nx+1:end);
            X22 = X(nx+1:end,nx+1:end);
            
% EXTRA ANALYSIS
            if analysis

% size of parameters                 
                param_norms(iter) = norm([cntrl_sys.A(:);cntrl_sys.B(:);cntrl_sys.C(:)],'fro');
                         
% convex paramters                
                [Kcvx,Lcvx,Mcvx,Xcvx,Ycvx,nu_e] = extended_Phi_map(cntrl_sys.A,cntrl_sys.B,cntrl_sys.C,X,sys);
                param_norms_cvx(iter) = norm(nu_e,'fro');  
                
% more specific parameter sizes 
                Kcvx_fro_norms(iter) = norm(Kcvx,'fro');
                Lcvx_fro_norms(iter) = norm(Lcvx,'fro');
                Mcvx_fro_norms(iter) = norm(Mcvx,'fro');
                Xcvx_fro_norms(iter) = norm(Xcvx,'fro');
                Ycvx_fro_norms(iter) = norm(Ycvx,'fro');
                
                Xinv = inv(X);
                Ucvx_2_norms(iter) = norm(Xinv(1:nx,nx+1:end));
                
                Vcvx_2_norms(iter) = norm(X(1:nx,nx+1:end));

                
% parameter error cvx space
                param_error_cvx(iter) = norm(nu_e - nu_e_opt,'fro');
                
% grad Phi
                dPhi = grad_extended_Phi_map_fd(cntrl_sys.A,cntrl_sys.B,cntrl_sys.C,X,sys);
                min_svd_dPhi(iter) = min(svd(dPhi));
                max_svd_dPhi(iter) = max(svd(dPhi));           
                
% Hessian
                Hess = kf_state_cost_ct_hessian_fd(ss2mat(cntrl_sys),sys);
                min_eig_hessian(iter) = min(eig(Hess));
                max_eig_hessian(iter) = max(eig(Hess));               
                
            end
                      
% Compute the cost and gradient            
            cs(iter) = c;
            grad_c_norms(iter) = norm(dc,2); 
            grad_c_fro_norms(iter) = norm(dc,'fro');
            
            min_sv_X12(iter) = min(svd(X12));
            min_eig_X22(iter) = min(eig(X22));
            error_eye_X22(iter) = norm(eye(nx) - X22,'fro');

            if reg_type == 0

                f = c;
                df = dc;

            elseif reg_type == 1

                [Rc,dRc] = cntrl_state_reg(cntrl,sys);

                f = c + lambda*Rc;
                df = dc + lambda*dRc;
                
                Rcs(iter) = Rc;
                dRc_norms(iter) = norm(dRc);

            elseif reg_type == 2

              
                [R12,dR12] = tr_inv_ct_kz(cntrl,sys);
                dR12 = [dR12; zeros(nx,nx+ny)];   
                
                f = c + lambda*(R12);
                df = dc + lambda*(dR12);
                
                R12s(iter) = R12;
                dR12_norms(iter) = norm(dR12);                

            elseif reg_type == 3

                [Rc,dRc] = cntrl_state_reg(cntrl,sys);

                [R12,dR12] = x_xinv_x_inv_reg(cntrl,sys);

                f = c + lambda*(lambda_c*Rc + R12);
                df = dc + lambda*(lambda_c*dRc + dR12);
                
                Rcs(iter) = Rc;
                dRc_norms(iter) = norm(dRc); 
                
                R12s(iter) = R12;
                dR12_norms(iter) = norm(dR12);   
                
            elseif reg_type == 4

                [Rc,dRc] = cntrl_state_reg(cntrl,sys);

                [R12,dR12] = x_xinv_x_inv_reg(cntrl,sys);

                f = Rc + R12;
                df = dRc + dR12;
                
                Rcs(iter) = Rc;
                dRc_norms(iter) = norm(dRc); 
                
                R12s(iter) = R12;
                dR12_norms(iter) = norm(dR12);  
                
% compute the optimal C_k   
                X12 = X(1:nx,nx+1:end);
                X22 = X(nx+1:end,nx+1:end);

                Copt = X12/X22;
                
% update the iterate                 
                cntrl_sys.C = Copt;
                cntrl = ss2mat(cntrl_sys);
                
% update the cost 
                [c, dc, ~] = kf_state_cost_ct(cntrl,sys);
                cs(iter) = c;
                grad_c_norms(iter) = norm(dc,2); 
                grad_c_fro_norms(iter) = norm(dc,'fro');
                
            elseif reg_type == 5

                [Ri,dRi] = state_cov_eye(cntrl,sys);

                [R12,dR12] = x_xinv_x_inv_reg(cntrl,sys);

                f = c + lambda*(lambda_c*Ri + R12);
                df = dc + lambda*(lambda_c*dRi + dR12);
                
                Ris(iter) = Ri;
                dRi_norms(iter) = norm(dRi); 
                
                R12s(iter) = R12;
                dR12_norms(iter) = norm(dR12);       
                
            elseif reg_type == 6
              
                [R12,dR12] = tr_inv_ct_kz(cntrl,sys);
                dR12 = [dR12; zeros(nx,nx+ny)];

                f = R12;
                df = dR12;
                                
                R12s(iter) = R12;
                dR12_norms(iter) = norm(dR12);  
                
% compute the optimal C_k   
                X12 = X(1:nx,nx+1:end);
                X22 = X(nx+1:end,nx+1:end);

                Copt = X12/X22;
                
% update the iterate                 
                cntrl_sys.C = Copt;
                cntrl = ss2mat(cntrl_sys);
                
% update the cost 
                [c, dc, ~] = kf_state_cost_ct(cntrl,sys);
                cs(iter) = c;
                grad_c_norms(iter) = norm(dc,2);
                grad_c_fro_norms(iter) = norm(dc,'fro');
                
            elseif reg_type == 7
                
                [R12,dR12] = tr_inv_ct_kz(cntrl,sys);
                dR12 = [dR12; zeros(nx,nx+ny)];   
                
%                 lambda = 1.05*norm(dc,'fro')/norm(dR12,'fro');
%                 lambda = 1e-6;
                
                f = c + lambda*(R12);
                df = dc + lambda*(dR12);
                
                R12s(iter) = R12;
                dR12_norms(iter) = norm(dR12);                 
                                                
            else 
                error('no such reg')
            end

            fs(iter) = f;
            grad_f_norms(iter) = norm(df,2);

%             fprintf('taking gradient step...\n')
            

% Perform a line search
            if linesearch == 1

                if reg_type == 7
                    
                    [t, cntrl, lsp] = simple_line_search(cntrl, f, df, -df, LS_gamma, LS_c, @(x) target(x,lambda));
                    stepsizes(iter) = t;                    
                    
                else
                
                    [t, cntrl, lsp] = simple_line_search(cntrl, f, df, -df, LS_gamma, LS_c, target);
                    stepsizes(iter) = t;
                    
                end
                

                if lsp %...  
                    
                    lsp_counter = lsp_counter + 1;
                    
                    if lsp_counter == 3
                    
                        if reg_type == 7
                            
% turn off regularization
                            if lambda == 0
                                
                                term_msg = 'line search problem';                    
                                run = 0;
                                
                            else
                                fprintf('turn off regularization, iter: %d\n',iter)
                                lambda = 0;
                                lsp_counter = 0;
                            end
                            
                        else
                                               
                            term_msg = 'line search problem';                    
                            run = 0;  
                            
                        end
                        
                    end
                    
                else
                    
                    lsp_counter = 0;
                        
                end
                
            elseif linesearch  == 2
                
                [t, cntrl, lsp] = exp_line_search(cntrl, f, df, -df, LS_gamma, LS_c, target);
                stepsizes(iter) = t;
                
                if lsp                   
                    term_msg = 'line search problem';
                    run = 0;                    
                end                
                

            else

                cntrl = cntrl - step_size*df;

            end
                       
            
% check for termination
            if iter == num_iters; term_msg = 'num iters'; run = 0; end

            if iter > 1
                
                if term_criterion == 1
                
                    if cs(iter-1) > cs(iter)
                        if cs(iter-1) - cs(iter) < obj_tol
                            term_msg = 'obj decrease tol';
                            run = 0;
                        end
                    end
                    
                elseif term_criterion == 2
                   
                    if grad_f_norms(iter) < grad_tol
                        term_msg = 'obj gradient tol';
                        run = 0;  
                    end
                                            
                end
                    
            end

            iter = iter + 1;
            
            
            if save_progress              
                    
                    if toc(last_report_time) > report_frequency_sec
                        
                        fprintf('Saving report at iter %d\n', iter)                        
                        
                        res.success = 1;
                        res.ops = ops;
                        res.cntrl_final = mat2ss(cntrl,nx);
                        res.term_msg = term_msg;
                        res.costs = cs;
                        res.reg_costs = fs;
                        res.grad_cost_norms = grad_c_norms;
                        res.grad_reg_cost_norms = grad_f_norms;
                        res.min_sv_X12 = min_sv_X12;
                        res.min_eig_X22 = min_eig_X22;
                        res.error_eye_X22 = error_eye_X22;
                        res.Rcs = Rcs;
                        res.dRc_norms = dRc_norms;
                        res.R12s = R12s;
                        res.dR12_norms = dR12_norms;
                        res.ctrbs = ctrbs;
                        res.obsvs = obsvs;
                        res.stepsizes = stepsizes;  
                        
                        if reg_type == 5
                            res.Ris = Ris;
                            res.dRi_norms = dRi_norms;
                        end
                        
                        if analysis

                            res.param_norms = param_norms;
                            res.param_norms_cvx = param_norms_cvx;
                            
                            res.param_error_cvx = param_error_cvx;
                            res.min_svd_dPhi = min_svd_dPhi;
                            res.max_svd_dPhi = max_svd_dPhi;
                            res.min_eig_hessian = min_eig_hessian;
                            res.max_eig_hessian = max_eig_hessian;

                        end                          
                            
                        figure
                        subplot(4,1,1)
                        loglog(cs)
                        title('cob cost')
                        
                        subplot(4,1,2)
                        loglog(grad_c_norms)
                        title('cob grad norm')                        
                        
                        subplot(4,1,3)
                        loglog(fs)
                        title('opt obj')  
                        
                        subplot(4,1,4)
                        loglog(grad_f_norms)
                        title('opt obj grad norm')    
                        
                        figure
                        loglog(error_eye_X22)
                        title('error |X22 - I|')
                        
                        save(fname, 'res')
                        pause(0.2)

                        last_report_time = tic;
                        
                    end
                                                
            end
            
            

        end
        
        fprintf('termination: %s\n', term_msg)
        
        res.success = 1;
        
        update_res();
                
                                        
    catch
        
        fprintf('some kind of failure!!!\n')
        
        res.success = 0;
        
        update_res();
       
        
        
    end


    function update_res()
        
        res.ops = ops;
        res.cntrl_final = mat2ss(cntrl,nx);
        res.term_msg = term_msg;
        res.costs = cs;
        res.reg_costs = fs;
        res.grad_cost_norms = grad_c_norms;
        res.grad_cost_fro_norms = grad_c_fro_norms;
        res.grad_reg_cost_norms = grad_f_norms;
        res.min_sv_X12 = min_sv_X12;
        res.min_eig_X22 = min_eig_X22;
        res.error_eye_X22 = error_eye_X22;
        res.Rcs = Rcs;
        res.dRc_norms = dRc_norms;
        res.R12s = R12s;
        res.dR12_norms = dR12_norms;
        res.ctrbs = ctrbs;
        res.obsvs = obsvs;
        res.stepsizes = stepsizes;
        
        if reg_type == 5
            res.Ris = Ris;
            res.dRi_norms = dRi_norms;
        end      
        
        if analysis

            res.param_norms = param_norms;
            res.param_norms_cvx = param_norms_cvx;

            res.param_error_cvx = param_error_cvx;
            res.min_svd_dPhi = min_svd_dPhi;
            res.max_svd_dPhi = max_svd_dPhi;
            res.min_eig_hessian = min_eig_hessian;
            res.max_eig_hessian = max_eig_hessian;
            
            res.Kcvx_norms = Kcvx_fro_norms;
            res.Lcvx_norms = Lcvx_fro_norms;
            res.Mcvx_norms = Mcvx_fro_norms;
            res.Xcvx_norms = Xcvx_fro_norms;
            res.Ycvx_norms = Ycvx_fro_norms;
            
            res.Ucvx_2_norms = Ucvx_2_norms;
            res.Vcvx_2_norms = Vcvx_2_norms;
            
            

        end         
       
    end



end
 