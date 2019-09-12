clear variables
%close all

set(groot, 'DefaultFigureRenderer', 'painters');
set(0, 'DefaultFigureRenderer', 'painters');

%
% Set up the overhead
%

[ a_par, rcm, data_c, fun_setup, fun_choose ] = setup_mmt_problem( );


if ~exist(a_par.plot_path, 'dir')
    mkdir(a_par.plot_path);
end
if ~exist(a_par.figure_path, 'dir')
    mkdir(a_par.figure_path);
end

a_par.save_intermediate_figs = true;


work_list = cell(a_par.max_runs, 1);

is_search_finished = false;


%
% Main Loop
%

[ r_par ] = fun_choose( a_par, rcm );

while (r_par.analysis_state == 0)
    r_par.analysis_state = 1;
    
    fprintf('Starting run %s.\n', r_par.title_tag);
    [h_c ] = fun_setup( data_c, a_par, r_par );
    
    w_base = W_Base();
    compute_basic_statistics( h_c, a_par, w_base );
    compute_prc_quantities( h_c, a_par, w_base );
    
    %compute_ad_quantities( pre, ind, NN, f_par, w_base );
    %compute_mi_quantities( pre, ind, NN, f_par, w_base );
    %compute_pca_quantities( pre, ind, NN, a_par, w_base );
    
    %
    % Plotting
    %
    
    if a_par.draw_scatter_plots
        tic;
        draw_scatter_plot( a_par, r_par, h_c);
        draw_prc_plot(a_par, r_par, w_base.w_prc);
        fprintf('Plotting done after %0.f seconds.\n', toc);
    end
    
    %
    % cleanup
    %
    
    work_list{r_par.run_number} = w_base;
    r_par.analysis_state = 2;
    
    [ r_par ] = fun_choose( a_par, rcm );
        
end

% switch project
%     case 'kolm'
%         draw_kolm_summaries( a_par, rcm, work_list);
%     case 'mmt'
%         draw_mmt_summaries( a_par, rcm, work_list);
%     case 'bimodal test'
%         draw_test_summaries( a_par, rcm, work_list);
%     case 'mutligauss test'
%         draw_test_summaries( a_par, rcm, work_list);
%     otherwise
%         warning('%s not recognized!', project);
% end
