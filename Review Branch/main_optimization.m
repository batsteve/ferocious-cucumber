clear variables
%close all

set(groot, 'DefaultFigureRenderer', 'painters');
set(0, 'DefaultFigureRenderer', 'painters');

opts = {'alpha_star', 'vus', 'ter', 'f_score'};

KK_num_files = [1];
KK_frac = {'half', 'full'};

num_runs = 1;

for k_opt = 1:length(opts)

    for k_files = 1:length(KK_num_files)

        x_set = cell(num_runs, 1);
        fval_set = cell(num_runs, 1);
        metric_set = cell(num_runs, 1);
        exitflag_set = cell(num_runs, 1);
        output_set = cell(num_runs, 1);
        trials_set = cell(num_runs, 1);

        for k_run = 1:num_runs
            fprintf('Starting %d-%d-%d.\n', k_opt, k_files, k_run);
            
            a_par = Analysis_Parameters();
            %addpath('../../MMT_NL/main_version');
            a_par.datapath = './MMT';
            a_par.gap_T_real = 0.015;
            a_par.analysis_name = sprintf('_t_%0.2f', a_par.gap_T_real);
            a_par.verbosity = 2;
            filename = sprintf('%s/sim_%d.dat', a_par.datapath, 1);
            fprintf('Loading file %s in order to read data parameters.\n', filename);

            load(filename, 'm_par', '-mat');

            a_par.drop_T_steps = 1000;
            a_par.gap_T_steps = ceil(a_par.gap_T_real./m_par.dt);
            a_par.gap_T_real = a_par.gap_T_steps*m_par.dt;
            a_par.analysis_name = sprintf('_t_%0.4f', a_par.gap_T_real);


            a_par.max_runs = 1000;
            a_par.max_time = 60*60*4;
            %a_par.max_time = 60*60*0.25;
            a_par.num_surrogate_points = 20;

            a_par.evaluation_mode = 'train';
            %a_par.nonlinear_constraint_type = 'penalty';
            a_par.nonlinear_constraint_type = 'none';
            a_par.fixed_A_threshold = 1.5;
            a_par.hypothesis_class_size = 2;

            num_files = KK_num_files(k_files);
            a_par.domain_frac = KK_frac{k_files};
            
            a_par.training_file_set = ...
                (num_files*k_run-(num_files - 1)):(num_files*k_run);
            a_par.testing_file_set = 23:30;

            a_par.opt_value = opts{k_opt};
            %a_par.opt_value = 'alpha_star';
            %a_par.opt_value = 'vus';
            %a_par.opt_value = 'ter';
            %a_par.opt_value = 'ber';
            %a_par.opt_value = 'f_score';


            a_par.save_intermediate_figs = false;
            a_par.draw_scatter_plots = false;
            a_par.plot_path = sprintf('./%s/', a_par.opt_value);
            a_par.figure_path = sprintf('./%s/', a_par.opt_value);
            if ~exist(a_par.plot_path, 'dir')
                mkdir(a_par.plot_path);
            end
            if ~exist(a_par.figure_path, 'dir')
                mkdir(a_par.figure_path);
            end


            fun_opt = @(coeff) fun_predictor_val(coeff, a_par);

            a_min = -1;
            a_max = 1;
            %L_min = log(0.5);
            %L_max = log(2000);
            L_min = 0;
            L_max = 1;

            n = a_par.hypothesis_class_size;
            lb = [a_min*ones(n - 1, 1); L_min*ones(n, 1)];
            ub = [a_max*ones(n - 1, 1); L_max*ones(n, 1)];

            a_par.search_algorithm = 'surrogateopt';

            switch a_par.search_algorithm
                case 'surrogateopt'

                    prob = struct;
                    prob.objective = fun_opt;
                    prob.lb = lb;
                    prob.ub = ub;
                    prob.solver = 'surrogateopt';
                    prob.options = optimoptions(@surrogateopt, 'MaxFunctionEvaluations', a_par.max_runs, ...
                        'Display', 'iter', 'PlotFcn','surrogateoptplot', ...
                        'MinSurrogatePoints', a_par.num_surrogate_points, ...
                        'MaxTime', a_par.max_time);
                    %prob.rngstate = []; % optional
                    %showproblem(prob);

                    [x,fval,exitflag,output,trials] = surrogateopt(prob);

                case 'pattern_search'

                    prob = struct;
                    prob.x0 = [0.25, 1/2*(L_min + L_max), 1/2*(L_min + L_max)];
                    prob.objective = fun_opt;
                    prob.Aineq = [];    % Use this to escape the symmetry?
                    prob.bineq = [];
                    prob.Aeq = [];
                    prob.Beq = [];
                    prob.lb = lb;
                    prob.ub = ub;
                    prob.nonlcon = [];
                    prob.solver = 'patternsearch';
                    prob.options = optimoptions(@patternsearch, 'MaxFunctionEvaluations', a_par.max_runs,...
                        'Cache', 'on');
                    %prob.rngstate = []; % optiona

                    [x,fval,exitflag,output] = patternsearch(prob);

                otherwise
                    warning('%s not recognized!\n', a_par.search_algorithm)
            end

            %a_par.evaluation_mode = 'test';
            %error_metrics = fun_predictor_val(x, a_par)

            fprintf('Run over!   Saving datas!\n\n\n');

            x_set{k_run} = x;
            fval_set{k_run} = fval;
            exitflag_set{k_run} = exitflag;
            output_set{k_run} = output;
            trials_set{k_run} = trials;
            %metric_set{k} = error_metrics;

            save(sprintf('./%s/partial_workspace_files_%d_%s_p_%d.mat', ...
                a_par.opt_value, KK_num_files(k_files), KK_frac{k_files}, k_run));
        end
        save(sprintf('./%s/workspace_files_long_%d_%s.mat', ...
            a_par.opt_value, KK_num_files(k_files), KK_frac{k_files}));
    end
end

