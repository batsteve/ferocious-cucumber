function [ a_par, rcm, data_c, fun_setup, fun_choose ] = setup_mmt_problem( )
%SETUP_MMT_PROBLEM Summary of this function goes here
%   Detailed explanation goes here

    a_par = Analysis_Parameters();
    a_par.plot_path = './plots/';
    a_par.figure_path = './figs/';
    
    
    %a_par.datapath = '../Kolmogorov Flow/SJ version/data/re_80';
    %addpath('../../MMT_NL/main_version');
    a_par.datapath = './MMT';
    a_par.num_files = 10;
    
    %a_par.gap_T_real = 0.0075;
    %a_par.gap_T_real = 1.5;
    a_par.gap_T_real = 0.015;
    %a_par.gap_T_real = 0.15;
    %a_par.gap_T_steps = ceil(a_par.gap_T_real*150/6000);
    %a_par.drop_T_real = 50;
    a_par.analysis_name = sprintf('_t_%0.2f', a_par.gap_T_real);
    a_par.verbosity = 2;
    
    %a_par.correspondance_dir = 'bijective';
    %a_par.correspondance_dir = 'fixed_predictor_box';
    %a_par.pre_extreme_mask_rule = 'none';
    %a_par.post_extreme_mask_rule = 'none';
    %a_par.mask_threshold = 2.5;
    
    %a_par.testing_gamut = 'single_L';
    a_par.testing_gamut = 'full_sweep';
    
    [ rcm ] = mmt_setup_control_matrix( a_par );
    %[ data_c ] = mmt_load_data( a_par );
    data_c = 0;
    
    %[ pre, ind, NN ] = mmt_setup_ind_pre( raw_sim, a_par, r_par );
    fun_setup = @(data_c, a_par, r_par) ...
        mmt_setup_ind_pre(data_c, a_par, r_par);
    fun_choose = @(a_par, runs_control_matrix) ...
        simple_choose_next_predictor( a_par, runs_control_matrix );
    
    
    
    filename = sprintf('%s/sim_%d.dat', a_par.datapath, 1);
    fprintf('Loading file %s (%d out of %d) in order to read data parameters.\n', filename, 1, a_par.num_files);
    
    load(filename, 'm_par', '-mat');
    
    a_par.drop_T_steps = 1000;
    a_par.gap_T_steps = ceil(a_par.gap_T_real./m_par.dt);
    a_par.gap_T_real = a_par.gap_T_steps*m_par.dt;
    %a_par.gap_T_steps = 10;
    a_par.analysis_name = sprintf('_t_%0.4f', a_par.gap_T_real);


end

