clear variables
%close all

set(groot, 'DefaultFigureRenderer', 'painters');
set(0, 'DefaultFigureRenderer', 'painters');

%opts = {'alpha_star', 'vus', 'ter', 'f_score'};
%opts = {'alpha_star', 'f_score'};

%KK_num_files = [10];
%KK_frac = {'full'};
num_runs = 1;
TT = 0:1:20;

for k_run = 1:num_runs
    for k_tau = 1:1:21
        
            fprintf('Starting %d-%d.\n', k_tau, k_run);
            
            a_par = Analysis_Parameters();
            a_par.model = 'kolm';
            a_par.datapath = './kolm_re_40';
            a_par.analysis_name = sprintf('_t_%0.2f', a_par.gap_T_real);
            a_par.verbosity = 2;
            
            %filename = sprintf('%s/sim_data/sim_%d.dat', a_par.datapath, 1);
            %fprintf('Loading file %s in order to read data parameters.\n', filename);
            %load(filename, 'm_par', '-mat');

            a_par.drop_T_steps = 0;
            a_par.gap_T_steps = TT(k_tau);


            a_par.max_runs = 4000;
            a_par.max_time = 60*60*4;
            a_par.num_surrogate_points = 64*2;

            a_par.evaluation_mode = 'train';
            a_par.nonlinear_constraint_type = 'none';
            a_par.fixed_A_threshold = 0.17;
            a_par.data_hist_bins = 128;

            a_par.domain_frac = 'full';
            
            a_par.data_length = 19900;
            a_par.training_file_set = 1;
            %a_par.testing_file_set = 1:1;

            a_par.opt_value = 'alpha_star';
            a_par.auc_integration_algorithm = 'trapezoid';
            a_par.bin_count_method = 'full_set';

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
            lb = [a_min*ones(64, 1)];
            ub = [a_max*ones(64, 1)];

            a_par.search_algorithm = 'surrogateopt';

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


            fprintf('Run over!   Saving datas!\n\n\n');

            mkdir(sprintf('./kt_%d', ...
                k_tau));
            save(sprintf('./workspace_t_%d_r_%d.mat', ...
                k_tau, k_run));
           
    end
end