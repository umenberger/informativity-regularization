%% dataset

% script uses data stored in 'results' variable

% this script will produce 4 figures:

% Figure 1 = does not appear in paper
% - individual traces from every trial plotted on the same figure

% Figure 2 = Figure 3(b) in paper
% - suboptimality as a function of iteration

% Figure 3 = Figure 4(a) in paper
% - informativity as a function of iteration

% Figure 4 = Figure 4(b) in paper
% - conditioning as a function of iteration

%% colors

colors = {[0, 0.4470, 0.7410],
          [0.8500, 0.3250, 0.0980],
          [0.4660, 0.6740, 0.1880],
          [0.9290, 0.6940, 0.1250],
          [0.4940, 0.1840, 0.5560],
          [0.4660, 0.6740, 0.1880],
          }; 

%% extract relevant data for plotting

% will plot all trials on the same axes

num_trials = length(results);

trials = 1:num_trials;

methods = {'gd','gdr','irpg'};

% populate legend
leg_names = {};

for m = 1:length(methods)       
    leg_names{end+1} = method_name(methods{m});       
end

% specify the things you want to plot in each figure, as a cell array
% things2plot = {[1]};
things2plot = {[1,2,9,10]};
      
% things to plot
% 1. suboptimality of cost (normalized)
% 2. norm of gradient of the cost
% 3. optimization objective (regularized)
% 4. norm of gradient of the objective (regularized)
% 5. value of covariance balancing reg
% 6. norm of gradient covariance balancing reg
% 7. value of trace inverse reg
% 8. norm of gradient of trace inverse reg
% 9. min svd of X12
% 10. min eig X22
% 11. | I - X22 |_F
% 12. stepsizes 
% 13. value of |X22 - I|^2
% 14. grad norm |X22 - I|^2

names_of_things = {'costs','grad_cost_norms','reg_costs','grad_reg_cost_norms',...
                   'Rcs', 'dRc_norms', 'R12s', 'dR12_norms', 'min_sv_X12',... 
                   'min_eig_X22', 'error_eye_X22', 'stepsizes',...
                   'Ris','dRi_norms'};

labels_of_things = {'suboptimality of cost','gradient of cost',...
                    'optimization objective (regularized)',...
                    'gradient of opt obj (regularized)',...
                    'cov bal cost', 'cov bal gradient', 'trace inv cost',...
                    'gradient of trace inv', 'min svd(X12)', 'min eig(X22)',...
                    '|I-X22|','stepsizes',...
                    'value |I-X22|^2','gradient of |I-X22|^2'};
               
num_figs = length(things2plot);
nsp = zeros(1,num_figs);
for i = 1:num_figs
    nsp(i) = length(things2plot{i});
end  

%
lw = 2.0;
transp = 0.3;
    
sp_counters = ones(1,num_figs);

figs = zeros(1,num_figs);
for i = 1:num_figs
    figs(i) = figure;
    p = get(gcf,'Position');     
    p(1) = p(1)*0.5;
    p(3) = p(3)*2;
    p(4) = p(4)*2.7;
    set(gcf, 'Position',  p)
end    

data_things = {};

for thing = 1:14

    for f = 1:num_figs

        if sum(things2plot{f} == thing) % check: plot the thing?

            figure(figs(f)) % open the figure

            subplot(nsp(f),1,sp_counters(f)) % open the subplot

            leg = {};
                        
            for m = 1:length(methods)
               
                loglog(nan, 'linewidth', lw, 'color', [colors{m} transp])                
                leg{end+1} = leg_names{m};
                
            end
            
            data_methods = {};
            
            for m = 1:length(methods) % attempt to plot for each method
                
                data_trials = [];
                
                for trial_index = trials

                    res = results{trial_index};    
                    f_opt = res.f_opt; 

                    if isfield(res.(methods{m}),names_of_things{thing}) % methods has thing  

                        data = real(res.(methods{m}).(names_of_things{thing}));

                        if thing == 1

                            data = (data - f_opt)/f_opt;

                        end
                                                
                        loglog(data, 'linewidth', lw, 'color', [colors{m} transp])
                        hold on   
                        
                        data_trials = [data_trials; data];

                    end

                end
                
                data_methods{m} = data_trials;

            end

            sp_counters(f) = sp_counters(f) + 1;
            legend(leg,'Location','SouthWest')
            legend boxoff
            grid on
%                 xlabel('iteration')
            if thing == 1
                if sum(strcmp(methods,'ctrb_X12_rebal'))

                    cntrl_up = update_c(res.ctrb_X12_rebal.cntrl_final,res.sys);
                    sub_up = kf_state_cost_ct(ss2mat(cntrl_up),res.sys);

                    title(['Trial ' num2str(trial_index) ', ' labels_of_things{thing}, ...
                           ', rebal after min_C: ' num2str((sub_up-f_opt)/f_opt,'%.4e')])
                else
                    title(['Trial ' num2str(trial_index) ', ' labels_of_things{thing}])                        
                end
            elseif thing == 3
                title([labels_of_things{thing}, ', opt cost: ' num2str(f_opt,'%.3f')]);                    
            else
                title(labels_of_things{thing})
            end                 

        end

        set(findall(gcf,'-property','FontSize'),'FontSize',15)

    end
    
    data_things{thing} = data_methods;

end




%% -----------------------------------------------------
%% Figure 3(b): suboptimality as a function of iteration
%% -----------------------------------------------------





%% cleaning for suboptimality 

data_methods = data_things{1};

% set a tolerance for optimality 
% if the final value (before a nan) is less than a certain tolerance,
% for plotting purposes it is considered optimal 
opt_tol = 1e-12;
opt_val = 1e-18;

data_methods_c = data_methods;

for m = 1:length(methods)
    
    for trial_index = 1:size(data_methods_c{m})
       
        nan_index = find(isnan(data_methods{m}(trial_index,:)),1);
        
        if data_methods{m}(trial_index,nan_index-1) < opt_tol
           
            fprintf('\t %d\n',trial_index)
            
            fin_val = opt_val;
            
        else
            
            fin_val = data_methods{m}(trial_index,nan_index-1);
            
        end
        
        data_methods_c{m}(trial_index,nan_index:end) = fin_val;
        
    end
    

end


%% percentiles (for suboptimality)

% pc_upr = 75; 
% pc_lwr = 25;

% uppers = {};
% lowers = {};

p95 = {};
p90 = {};
p75 = {};
p25 = {};
p10 = {};
p5 = {};

uplyers = {};
downlyers = {}; 

for m = 1:length(methods)
   
    p95{m} = prctile(data_methods_c{m},95);
    p90{m} = prctile(data_methods_c{m},90);
    p75{m} = prctile(data_methods_c{m},75);
    p25{m} = prctile(data_methods_c{m},25);
    p10{m} = prctile(data_methods_c{m},10);
    p5{m} = prctile(data_methods_c{m},5);
        
%   'outliers'
    uplyers{m} = find((sum(data_methods_c{m} > p95{m},2)) >  0.01*length(p95{m}));
    downlyers{m} = find((sum(data_methods_c{m} < p5{m},2)) > 0.01*length(p95{m}));
        
end

uplyers
downlyers

%% plot suboptimality

figure

% for the legend
legend_names = {};

for m = 1:length(methods) 
    loglog(nan,'linewidth', 4, 'color', [colors{m} 0.4])
    hold on
    legend_names{end+1} = method_name_paper(methods{m});       
end

% plot these percentiles 
upr = p75;
lwr = p25;

Upr = p90;
Lwr = p10;

for m = 1:length(methods)

    if 1

    num = length(Upr{m});
    xs = [1:99, round(logspace(2,log10(num),1e3))];
    
    if sum(isnan(Upr{m}(xs)))
       
        upr_nnan = find(isnan(Upr{m}(xs)),1) - 1;
        
        tmp = find(isnan(Upr{m}),1) - 1;
        
        xsu = [xs(1:upr_nnan) tmp];
        
        
    else
        
        xsu = xs;
        
    end
    
    if sum(isnan(Lwr{m}(xs)))
       
        lwr_nnan = find(isnan(Lwr{m}(xs)),1) - 1;
        tmp = find(isnan(Lwr{m}),1) - 1;
        
        xsl = fliplr([xs(1:lwr_nnan) tmp]); 
        
    else
        
        xsl = fliplr(xs);
        
    end    
    
    upd = Upr{m}(xsu);
    
    lwd = Lwr{m}(xsl);
    
    patch([xsu xsl], ...
          [upd lwd], ...
          [1 1 1], ...
          'EdgeColor',[1 1 1], ...
          'EdgeAlpha',0, ...
          'FaceColor', colors{m}, ...
          'FaceAlpha', 0.1) 
      
      second_alpha = 0.25;
      
    else
        
        second_alpha = 0.2;

    end

% -------------------------------------------------------------------------

    num = length(upr{m});
    xs = [1:99, round(logspace(2,log10(num),1e3))];
    
    if sum(isnan(upr{m}(xs)))
       
        upr_nnan = find(isnan(upr{m}(xs)),1) - 1;
        
        tmp = find(isnan(upr{m}),1) - 1;
        
        xsu = [xs(1:upr_nnan) tmp];
        
        
    else
        
        xsu = xs;
        
    end
    
    if sum(isnan(lwr{m}(xs)))
       
        lwr_nnan = find(isnan(lwr{m}(xs)),1) - 1;
        tmp = find(isnan(lwr{m}),1) - 1;
        
        xsl = fliplr([xs(1:lwr_nnan) tmp]); 
        
    else
        
        xsl = fliplr(xs);
        
    end    
    
    upd = upr{m}(xsu);
    
    lwd = lwr{m}(xsl);
    
    patch([xsu xsl], ...
          [upd lwd], ...
          [1 1 1], ...
          'EdgeColor',[1 1 1], ...
          'EdgeAlpha',0, ...
          'FaceColor', colors{m}, ...
          'FaceAlpha', second_alpha)
      
                
% plot the outliers 
%     outliers = unique([uplyers{m}; downlyers{m}])';
    outliers = unique([downlyers{m}])';
                
end

legend(legend_names,'interpreter','latex')
legend boxoff

xlabel('Iteration','interpreter','latex')
ylabel('Suboptimality $\frac{\mathcal{L}_\mathtt{OE} - \mathcal{L}_\mathtt{OE}^\star}{\mathcal{L}_\mathtt{OE}^\star}$','interpreter','latex')

YLim = get(gca,'YLim');
set(gca,'YLim',[1e-15, 1e3])

grid on

set(gca,'XMinorGrid','off')
% set(gca,'YMinorGrid','off')

set(findall(gcf,'-property','FontSize'),'FontSize',17)
set(gca,'TickLabelInterpreter','latex')




%% -----------------------------------------------------
%% Figure 4(a): informativity as a function of iteration
%% -----------------------------------------------------



%% cleaning for "informativity" 

data_methods = data_things{9}; % min(svd(Sigma_12))

% data "cleaning"
% e.g. once nan, extend with previous value
data_methods_c = data_methods;

for m = 1:length(methods)
    
    
    for trial_index = 1:size(data_methods_c{m})
       
        nan_index = find(isnan(data_methods{m}(trial_index,:)),1);
        
        fin_val = data_methods{m}(trial_index,nan_index-1);
                
        data_methods_c{m}(trial_index,nan_index:end) = fin_val;
        
    end
    

end


%% percentiles

p95 = {};
p90 = {};
p75 = {};
p25 = {};
p10 = {};
p5 = {};

uplyers = {};
downlyers = {}; 

for m = 1:length(methods)
   
    p95{m} = prctile(data_methods_c{m},95);
    p90{m} = prctile(data_methods_c{m},90);
    p75{m} = prctile(data_methods_c{m},75);
    p25{m} = prctile(data_methods_c{m},25);
    p10{m} = prctile(data_methods_c{m},10);
    p5{m} = prctile(data_methods_c{m},5);
        
    uplyers{m} = find((sum(data_methods_c{m} > p95{m},2)) >  0.01*length(p95{m}));
    downlyers{m} = find((sum(data_methods_c{m} < p5{m},2)) > 0.01*length(p95{m}));
        
end

uplyers
downlyers

%% plot upper and lower (shaded) 

figure

% for the legend
legend_names = {};

method_selection = 1:length(methods); % all 

% for m = 1:length(methods) 
for m = method_selection
    loglog(nan,'linewidth', 4, 'color', [colors{m} 0.4])
    hold on
    legend_names{end+1} = method_name_paper(methods{m});       
end

upr = p75;
lwr = p25;

Upr = p90;
Lwr = p10;

for m = method_selection

    if 1

    num = length(Upr{m});
    xs = [1:99, round(logspace(2,log10(num),1e3))];
    
    if sum(isnan(Upr{m}(xs)))
       
        upr_nnan = find(isnan(Upr{m}(xs)),1) - 1;
        
        tmp = find(isnan(Upr{m}),1) - 1;
        
        xsu = [xs(1:upr_nnan) tmp];
        
        
    else
        
        xsu = xs;
        
    end
    
    if sum(isnan(Lwr{m}(xs)))
       
        lwr_nnan = find(isnan(Lwr{m}(xs)),1) - 1;
        tmp = find(isnan(Lwr{m}),1) - 1;
        
        xsl = fliplr([xs(1:lwr_nnan) tmp]); 
        
    else
        
        xsl = fliplr(xs);
        
    end    
    
    upd = Upr{m}(xsu);
    
    lwd = Lwr{m}(xsl);
    
    patch([xsu xsl], ...
          [upd lwd], ...
          [1 1 1], ...
          'EdgeColor',[1 1 1 ], ...
          'EdgeAlpha', 0, ...
          'FaceColor', colors{m}, ...
          'FaceAlpha', 0.1) 
      
      second_alpha = 0.3;
      
    else
        
        second_alpha = 0.3;

    end

% -------------------------------------------------------------------------

    num = length(upr{m});
    xs = [1:99, round(logspace(2,log10(num),1e3))];
    
    if sum(isnan(upr{m}(xs)))
       
        upr_nnan = find(isnan(upr{m}(xs)),1) - 1;
        
        tmp = find(isnan(upr{m}),1) - 1;
        
        xsu = [xs(1:upr_nnan) tmp];
        
        
    else
        
        xsu = xs;
        
    end
    
    if sum(isnan(lwr{m}(xs)))
       
        lwr_nnan = find(isnan(lwr{m}(xs)),1) - 1;
        tmp = find(isnan(lwr{m}),1) - 1;
        
        xsl = fliplr([xs(1:lwr_nnan) tmp]); 
        
    else
        
        xsl = fliplr(xs);
        
    end    
    
    upd = upr{m}(xsu);
    
    lwd = lwr{m}(xsl);
    
    patch([xsu xsl], ...
          [upd lwd], ...
          [1 1 1], ...
          'EdgeColor',[1 1 1], ...
          'EdgeAlpha', 0, ...          
          'FaceColor', colors{m}, ...
          'FaceAlpha', second_alpha)
      
                
% plot the outliers 
    outliers = unique([uplyers{m}; downlyers{m}])';
%     outliers = unique([downlyers{m}])';

    
    for trial_index = outliers
        
%         loglog(data_methods{m}(trial_index,:), 'linewidth', lw, 'color', [colors{m} 0.1])
        
    end
                
end



legend(legend_names,'interpreter','latex','Location','Southwest')
legend boxoff

xlabel('Iteration','interpreter','latex')
ylabel('Informativity, $\sigma_{\min}(\mathbf{\Sigma}_{12})$','interpreter','latex')
% set(gca,'YTick',[1e-10,1e-5,1])

% YLim = get(gca,'YLim');
% set(gca,'YLim',[1e-15, 1e3])
grid on
set(gca,'XMinorGrid','off')
set(gca,'YMinorGrid','off')

set(findall(gcf,'-property','FontSize'),'FontSize',17)
set(gca,'TickLabelInterpreter','latex')

% set(gca,'YTick',[1e-10,1e-8,1e-6,1e-4,1e-2])





%% -----------------------------------------------------
%% Figure 4(b): conditioning as a function of iteration
%% -----------------------------------------------------



%% cleaning for "informativity" 

data_methods = data_things{10}; % min(eig(Sigma_22))

% data "cleaning"
% e.g. once nan, extend with previous value
data_methods_c = data_methods;

for m = 1:length(methods)
    
    
    for trial_index = 1:size(data_methods_c{m})
       
        nan_index = find(isnan(data_methods{m}(trial_index,:)),1);
        
        fin_val = data_methods{m}(trial_index,nan_index-1);
                
        data_methods_c{m}(trial_index,nan_index:end) = fin_val;
        
    end
    

end



%% percentiles

% if a trial is an outlier at any stage, plot the entire outlier?

p95 = {};
p90 = {};
p75 = {};
p25 = {};
p10 = {};
p5 = {};

uplyers = {};
downlyers = {}; 

for m = 1:length(methods)
   
    p95{m} = prctile(data_methods_c{m},95);
    p90{m} = prctile(data_methods_c{m},90);
    p75{m} = prctile(data_methods_c{m},75);
    p25{m} = prctile(data_methods_c{m},25);
    p10{m} = prctile(data_methods_c{m},10);
    p5{m} = prctile(data_methods_c{m},5);
        
    uplyers{m} = find((sum(data_methods_c{m} > p95{m},2)) >  0.01*length(p95{m}));
    downlyers{m} = find((sum(data_methods_c{m} < p5{m},2)) > 0.01*length(p95{m}));
        
end

uplyers
downlyers

%% plot upper and lower (shaded) 

figure

% for the legend
legend_names = {};

method_selection = 1; % only unregularized gradient descent

% for m = 1:length(methods) 
for m = method_selection
    loglog(nan,'linewidth', 4, 'color', [colors{m} 0.4])
    hold on
    legend_names{end+1} = method_name_paper(methods{m});       
end

upr = p75;
lwr = p25;

Upr = p90;
Lwr = p10;

for m = method_selection


    if 1

    num = length(Upr{m});
    xs = [1:99, round(logspace(2,log10(num),1e3))];
    
    if sum(isnan(Upr{m}(xs)))
       
        upr_nnan = find(isnan(Upr{m}(xs)),1) - 1;
        
        tmp = find(isnan(Upr{m}),1) - 1;
        
        xsu = [xs(1:upr_nnan) tmp];
        
        
    else
        
        xsu = xs;
        
    end
    
    if sum(isnan(Lwr{m}(xs)))
       
        lwr_nnan = find(isnan(Lwr{m}(xs)),1) - 1;
        tmp = find(isnan(Lwr{m}),1) - 1;
        
        xsl = fliplr([xs(1:lwr_nnan) tmp]); 
        
    else
        
        xsl = fliplr(xs);
        
    end    
    
    upd = Upr{m}(xsu);
    
    lwd = Lwr{m}(xsl);
    
    patch([xsu xsl], ...
          [upd lwd], ...
          [1 1 1], ...
          'EdgeColor',[1 1 1 ], ...
          'EdgeAlpha', 0, ...
          'FaceColor', colors{m}, ...
          'FaceAlpha', 0.1) 
      
      second_alpha = 0.3;
      
    else
        
        second_alpha = 0.3;

    end

% -------------------------------------------------------------------------

    num = length(upr{m});
    xs = [1:99, round(logspace(2,log10(num),1e3))];
    
    if sum(isnan(upr{m}(xs)))
       
        upr_nnan = find(isnan(upr{m}(xs)),1) - 1;
        
        tmp = find(isnan(upr{m}),1) - 1;
        
        xsu = [xs(1:upr_nnan) tmp];
        
        
    else
        
        xsu = xs;
        
    end
    
    if sum(isnan(lwr{m}(xs)))
       
        lwr_nnan = find(isnan(lwr{m}(xs)),1) - 1;
        tmp = find(isnan(lwr{m}),1) - 1;
        
        xsl = fliplr([xs(1:lwr_nnan) tmp]); 
        
    else
        
        xsl = fliplr(xs);
        
    end    
    
    upd = upr{m}(xsu);
    
    lwd = lwr{m}(xsl);
    
    patch([xsu xsl], ...
          [upd lwd], ...
          [1 1 1], ...
          'EdgeColor',[1 1 1], ...
          'EdgeAlpha', 0, ...          
          'FaceColor', colors{m}, ...
          'FaceAlpha', second_alpha)
      
                
% plot the outliers 
    outliers = unique([uplyers{m}; downlyers{m}])';
%     outliers = unique([downlyers{m}])';

    
    for trial_index = outliers
        
%         loglog(data_methods{m}(trial_index,:), 'linewidth', lw, 'color', [colors{m} 0.1])
        
    end
                
end


legend(legend_names,'interpreter','latex','Location','Southwest')
legend boxoff

xlabel('Iteration','interpreter','latex')
ylabel('Conditioning, $\lambda_{\min}(\mathbf{\Sigma}_{22})$','interpreter','latex')
set(gca,'YTick',[1e-8,1e-6,1e-4,1e-2,1])

% YLim = get(gca,'YLim');
% set(gca,'YLim',[1e-15, 1e3])
grid on
set(gca,'XMinorGrid','off')
set(gca,'YMinorGrid','off')

set(findall(gcf,'-property','FontSize'),'FontSize',17)
set(gca,'TickLabelInterpreter','latex')

% set(gca,'YTick',[1e-10,1e-8,1e-6,1e-4,1e-2])












