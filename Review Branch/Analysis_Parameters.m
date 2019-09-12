classdef Analysis_Parameters < handle
    %ANALYSIS_PARAMTERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % general model stuff
        model;
        
        % general algorithm stuff
        eps;
        
        % optimization control
        max_runs;
        data_path;
        %num_files;
        testing_file_set;
        training_file_set
        domain_frac;
        
        % more optimization control
        opt_value;
        cur_run;
        parameter_scale;
        parameter_unroll_coding;
        search_algorithm;
        hypothesis_class_size;
        evaluation_mode;
        nonlinear_constraint_type;
        num_surrogate_points;
        max_time;
        
        % file io and histogram creation
        datapath;
        savepath;
        %n_files_to_load;
        %dt_1;
        %dt_2;
        %drop_count;
        correspondance_dir;
        normalize_raw_data;
        gap_T_real;
        gap_T_steps;
        drop_T_real;
        drop_T_steps;
        predictor_data_choice;
        data_hist_bins;
        pre_extreme_mask_rule;
        post_extreme_mask_rule;
        mask_threshold;
        use_identical_histograms;
        save_workspace;
        data_length;
        
        % mmt
        kernel_shape;
        n_x;
        testing_gamut;
        
        % kolm
        fourier_set;
        fourier_variable;
        
        % test
        num_counts;
        
        % prc
        threshold_distribution;
        n_a_thresholds;
        n_b_thresholds;
        auc_integration_algorithm;
        auc_integration_points;
        auc_plot_points;
        bin_count_method;
        
        % binary
        fixed_A_threshold;
        F_score_beta;
        
        % troubleshooting
        verbosity;
        
        % plotting
        draw_scatter_plots;
        save_intermediate_figs;
        save_summary_figs
        plot_path;
        figure_path;
        analysis_name;
        renderer;
        next_scatter_plot_num;
        next_prc_plot_num;
        next_misc_plot_num;
        iterate_plots;
    end
    
    methods
        function [ a_par ] = Analysis_Parameters()
            a_par.model = 'mmt';
            %a_par.model = 'kolm';
            
            a_par.eps = 1e-6;
            
            a_par.max_runs = 50;
            a_par.data_path = '';
            %a_par.num_files = 10;
            a_par.testing_file_set = 1:10;
            a_par.training_file_set = 1:10;
            a_par.domain_frac = 'full';
            %a_par.domain_frac = 'half';
            
            a_par.opt_value = 'alpha_star';
            a_par.cur_run = 0;
            a_par.parameter_scale = 'log';
            %a_par.parameter_unroll_coding = 'independent';
            a_par.parameter_unroll_coding = 'stacked';
            a_par.search_algorithm = 'surrogateopt';
            a_par.hypothesis_class_size = 1;
            a_par.evaluation_mode = 'train';
            a_par.nonlinear_constraint_type = 'none';
            %a_par.nonlinear_constraint_type = 'penalty';
            a_par.num_surrogate_points = 20;
            a_par.max_time = 60*60*24;
            
            %a_par.n_files_to_load = 10;
            %a_par.dt_1 = 0;
            %a_par.dt_2 = 0;
            %a_par.drop_count = 0;
            a_par.correspondance_dir = 'bijective';
            a_par.datapath = '';
            a_par.savepath = '';
            a_par.normalize_raw_data = true;
            a_par.gap_T_real = 0;
            a_par.gap_T_steps = 0;
            a_par.drop_T_real = 0;
            a_par.drop_T_steps = 0;
            a_par.predictor_data_choice = '';
            a_par.data_hist_bins = 1024;
            a_par.pre_extreme_mask_rule = 'none';
            a_par.post_extreme_mask_rule = 'none';
            a_par.mask_threshold = 0;
            a_par.use_identical_histograms = true;
            a_par.save_workspace = false;
            a_par.data_length = 6000;
            
            a_par.kernel_shape = 'exp_squared';
            a_par.n_x = 8192;
            a_par.testing_gamut = 'full_sweep';
            
            %a_par.fourier_set = 'full';
            %a_par.fourier_set = 'half';
            a_par.fourier_set = 'quadrant';
            a_par.fourier_variable = 'u_12';
            
            a_par.num_counts = 1e8;
            
            a_par.n_a_thresholds = 256;
            a_par.n_b_thresholds = 256;
            a_par.threshold_distribution = 'linear';
            a_par.auc_integration_algorithm = 'scattered interpolant';
            %a_par.auc_integration_algorithm = 'trapezoid';
            a_par.auc_integration_points = 513;
            a_par.auc_plot_points = 101;
            a_par.bin_count_method = 'full_set';
            %a_par.bin_count_method = 'one-by-one';
            
            a_par.fixed_A_threshold = 1.8;
            a_par.F_score_beta = 1;
            
            a_par.verbosity = 1;
            
            a_par.draw_scatter_plots = true;
            a_par.save_intermediate_figs = true;
            a_par.save_summary_figs = true;
            a_par.plot_path = 'plots/';
            a_par.figure_path = 'figs/';
            a_par.analysis_name = '';
            a_par.renderer = 'painters';
            a_par.next_scatter_plot_num = 21;
            a_par.next_prc_plot_num = 41;
            a_par.next_misc_plot_num = 61;
            a_par.iterate_plots = false;
        end
    end
    
end

