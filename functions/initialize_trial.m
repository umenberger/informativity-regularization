function [sys,cntrl_opt,cntrl_init,f_opt] = initialize_trial(nx,nu,ny)

% Function outputs:
% sys - true linear system
% cntrl_opt - optimal filter (poorly named as "cntrl")
% cntrl_init - initial filter (to initialize policy gradient)
% f_opt - value/cost of optimal filter


    % Rejection sampling for true system
    min_min_true_obsv = 1e-4;
    max_min_true_obsv = 1e-2;
    rejection_meas_noise = 0;
    rejection_state_disturbance = 1;
    max_opt_cost = 1e3;
    rejection_obsv = 1;

    % Rejection sampling for initial system
    max_initial_subopt_factor = 100;

    min_initial_min_svd_X12 = 1e-5;
    max_initial_min_svd_X12 = 1e-3;

    min_initial_min_eig_X22 = 1e-3;
    max_initial_min_eig_X22 = 1;

    initial_distance_cvx = 1e3;
    initial_attempts = 1e3;


    suitable_trial = 0;
    
    fprintf('Rejection sampling for suitable initialization...\n')

    while ~suitable_trial

        suitable = 0;

        while ~suitable

            try

    % rejection sample observability 
                if rejection_obsv

                    unsuitable = 1;

                    while unsuitable

                        tmp_sys = rss(nx,ny,nu);
                        gram_o = gram(tmp_sys,'o');                        

                        if (min_min_true_obsv < min(eig(gram_o))) && (min(eig(gram_o)) < max_min_true_obsv)

                            unsuitable = 0;

                        end

                    end

                end

    % simply set measurement noise variance to unity 
                V = eye(ny);
    %                 V = randn(ny); V = V'*V; 

    % rejection sample disturbance variance
                if rejection_state_disturbance

                    unsuitable = 1;

                    while unsuitable

                        W = randn(nx); W = W'*W;

                        if max(eig(W)) < 5
                            unsuitable = 0;
                        end

                    end
                end

                A = tmp_sys.A; 
                B = tmp_sys.B; 
                C = tmp_sys.C;

                sys.A = A;
                sys.B = zeros(nx,nx); % important for the regularizers! 
                sys.C = C;
                sys.W = W;
                sys.V = V;  

                [cntrl_opt,P,L] = kf_via_care(A,C,W,V);
                cntrl_opt.C = eye(nx);
                cntrl_opt.D = zeros(nx,ny);

                f_opt = kf_state_cost_ct(ss2mat(cntrl_opt),sys); % optimal cost

    % rejection sample based on cost                 

                if f_opt < max_opt_cost

                    suitable = 1;

                end

            catch

%                 fprintf('problem during the search for true system\n')

            end

        end

%         fprintf('Found suitable true system and optimal filter\n')

        fltr_opt = ss2mat(cntrl_opt);

    % search for initial system
    % -------------------------------------------------------------------------        

        suitable = 0;

        init_iter = 0;

        try 

            while ~suitable 

                if init_iter == initial_attempts
                    break
                end

                unsuitable = 1;

                while unsuitable

                    cntrl_init = fltr_opt + 10*randn(size(fltr_opt));

                    [c,~,X] = kf_state_cost_ct(cntrl_init,sys);

                    X12 = X(1:nx,nx+1:end);
                    X22 = X(nx+1:end,nx+1:end);

                    min_svd_X12 = min(svd(X12));
                    min_eig_X22 = min(eig(X22));

                    if  (max_initial_min_svd_X12 > min_svd_X12) && ...
                        (min_svd_X12 > min_initial_min_svd_X12) && ...
                        (max_initial_min_eig_X22 > min_eig_X22) && ...
                        (min_eig_X22 > min_initial_min_eig_X22) 

                        unsuitable = 0;

                    end

                end


    % check for suboptimality    

                if kf_state_cost_ct(cntrl_init,sys) < max_initial_subopt_factor*f_opt

                        suitable = 1;
                        suitable_trial = 1;                      

                end 

                init_iter = init_iter + 1;

            end

        catch

%             fprintf('problem during the search for initial filter\n')

        end

        cntrl_init = mat2ss(cntrl_init,nx);        

    end

    fprintf('Found suitable initialization\n')


end