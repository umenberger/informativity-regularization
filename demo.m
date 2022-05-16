%% instructions

% The following values have been modified compared to the values used for
% the experiments in the paper. This has been done to allow the script to
% run in a "reasonable" amount of time.

% To run the script quickly use these values:
num_trials = 5; % number of independent trials
num_iters = 100; % number of iterations of gradient descent

% To REPRODUCE THE FIGURES in the paper use these values
% num_trials = 50; % value used in the paper's experiments
% num_iters = 1e5; % value used in the paper's experiments

%% setup

% methods to run
run_gd = 1;    % unregularized gradient descent
run_gdr = 1;   % gradient descent with reconditioning step
run_irpg = 1;  % information-regularized policy gradient (proposed method)

% system dimensions
nx = 2;
nu = 1;
ny = 1;

% results of all experiments will be stored in this variable
results = {};

% save output of these computations?
save_experiment = 0; % toggle to 1 to save

if save_experiment
   
    dt = datestr(now,'yyyy_mm_dd_HH_MM_SS');
    fname = [dt '_experiments'];    
    save(fname, 'results');
    
end

%% run trials

for trial_index = 1:num_trials
   
    fprintf('\n\nTrial %d\n', trial_index)
    fprintf('------------------------------------------------\n')
    
% initialize search: 
    [sys,cntrl_opt,cntrl_init,f_opt] = initialize_trial(nx,nu,ny);  
    
% store results  
    clear trial_result
    trial_result.cntrl_opt = cntrl_opt;
    trial_result.sys = sys;
    trial_result.f_opt = f_opt;
    trial_result.cntrl_init = cntrl_init;
       
% run the methods
% ------------------------------------------------------------------------- 
    
    plot_after_each_method = 1;
    findex = figure();
    leg = {};
    lw = 1.5;


    if run_gd
        
        fprintf('\nMethod: unregularized gradient descent, started: %s\n',datestr(now,'HH:MM:SS'))
        fprintf('---------------------------\n')
        
        clear ops        
        ops.reg_type = 0;
        ops.lambda = 0;
        ops.save_progress = 0;
        ops.linesearch = 1;
        ops.rebalance_X22 = 0;
        ops.termination_criterion = 'gradient';
        ops.grad_tol = 1e-8;
        ops.num_iters = num_iters;
        ops.analysis = 1;
       
        res_reg = gd_filtering_state_cost(cntrl_init, sys, ops);
        
        if plot_after_each_method && res_reg.success

            figure(findex)
            loglog(1:length(res_reg.costs),(res_reg.costs-f_opt)/f_opt, 'linewidth', lw)
            hold on

            pause(0.1)

            leg{end+1} = 'gd';
        
        end
        
        trial_result.gd = res_reg;
          
    end 
    
    
    if run_gdr
        
        fprintf('\nMethod: gradient descent with reconditioning, started: %s\n',datestr(now,'HH:MM:SS'))
        fprintf('---------------------------\n')        
        
        clear ops        
        ops.reg_type = 0;
        ops.lambda = 0;
        ops.save_progress = 0;
        ops.linesearch = 1;
        ops.rebalance_X22 = 1;
        ops.termination_criterion = 'gradient';
        ops.grad_tol = 1e-8;
        ops.num_iters = num_iters;
        ops.analysis = 1;
       
        res_reg = gd_filtering_state_cost(cntrl_init, sys, ops);
        
        if res_reg.success

            figure(findex)
            loglog(1:length(res_reg.costs),(res_reg.costs-f_opt)/f_opt, 'linewidth', lw)
            hold on

            pause(0.1)

            leg{end+1} = 'gdr';
        
        end
        
        trial_result.gdr = res_reg;
          
    end       

    
    if run_irpg
        
        fprintf('\nMethod: informativity-regularized policy gradient, started: %s\n',datestr(now,'HH:MM:SS'))
        fprintf('---------------------------\n')        
        
        clear ops
        ops.reg_type = 7;
        ops.lambda = 1e-4;
        ops.save_progress = 0;
        ops.linesearch = 1;
        ops.rebalance_X22 = 1;
        ops.termination_criterion = 'gradient';  
        ops.grad_tol = 1e-8;
        ops.num_iters = num_iters;
        ops.analysis = 1;
       
        res_reg = gd_filtering_state_cost(cntrl_init, sys, ops);
        
        if res_reg.success

            figure(findex)
            loglog(1:length(res_reg.costs),(res_reg.costs-f_opt)/f_opt, 'linewidth', lw)
            hold on
            pause(0.1)

            leg{end+1} = 'irpg';
        
        end
        
        trial_result.irpg = res_reg;
          
    end     
    
    legend(leg)
    xlabel('Iteration')
    ylabel('Output Error cost')


% store results 
% ------------------------------------------------------------------------- 

    if save_experiment 
        
        load(fname);

        results{trial_index} = trial_result;

        clear trial_result

        save(fname,'results');

        pause(1)

        clear results        
                
    else
        
        results{trial_index} = trial_result;
        
    end

end


%% plot the results 

if save_experiment 

    load(fname);

end

convergence_plots






