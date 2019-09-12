function [ val ] = fun_predictor_val( coeffs, a_par )
%FUN_PREDICTOR_VAL Summary of this function goes here
%   Detailed explanation goes here

    tic;
    
    r_par = Run_Parameters();
    a_par.cur_run = a_par.cur_run + 1;
    r_par.run_number = a_par.cur_run;
    r_par.coeffs = coeffs;
    r_par.file_tag = sprintf('_n=%d', r_par.run_number);
    r_par.title_tag = sprintf('(n=%d)', r_par.run_number);
    
    
    
    fprintf('Starting run %s:\n', r_par.title_tag);
    
    switch a_par.model
        case 'mmt'
            [ AA, LL ] = parametersToWavelets( a_par, coeffs );
            r_par.pred_fun = Predictor_Container(a_par, coeffs);
            for k = 1:a_par.hypothesis_class_size
                fprintf('a%d=%0.2f, L%d=%0.2f.\n', k, AA(k), k, LL(k));
            end
            
        case 'kolm'
            total_amp = sum(coeffs(:).^2);
            if (total_amp < a_par.eps)
                fprintf('This is the null predictor!\n');
                fprintf('Emergency workaround for PRC quantities invoked!\n');
                r_par.dangerous_prc = true;
            else
                % add a cool way to show major coeffs?
            end
    end

    
    

    switch a_par.nonlinear_constraint_type 
        % use this to enforce certain kinds of constraints on the
        % hypothesis space without expensive evaluations
        case 'penalty'
            
            % keep coefficients reasonable
            if (real(AA(1)) < 0) || abs(imag(AA(1)) > 0)
                val = 1;
                fprintf('Residual coefficient is %0.2f+%0.2fi.  Aborting evaluation.\n', real(AA(1)), imag(AA(1)));
                return
            end
            
            % keep wavelets ordered
            for k = 2:num_wavelets
                if (LL(k-1) > (LL(k) + 0.01))
                    val = 1;
                    fprintf('Wavelets disordered.  Aborting evaluation\n');
                    dangerous_prc = true;
                end
            end
            
        case 'none'
            
        otherwise
    end

    h_c = Histogram_Container(a_par, r_par);
    
    w_base = W_Base();
    %compute_basic_statistics( h_c, a_par, w_base );
    
    calc_bc = false;
    calc_prc = false;
    
    switch a_par.evaluation_mode
        case 'train'
            switch a_par.opt_value
                case 'alpha_star', calc_prc = true;
                case 'vus',    calc_prc = true;
                case 'ter',    calc_bc = true;
                case 'ber',    calc_bc = true;
                case 'f_score',calc_bc = true;

                otherwise
                    warning('%s not recognized!', a_par.opt_value);
            end
            
        case 'test'
            calc_bc = true;
            calc_prc = true;
            
        otherwise
            warning('%s not recognized!\n', a_par.evaluation_mode);
    end
    
    if calc_bc
        compute_bc_quantities( h_c, a_par, w_base );
    end
    if calc_prc
        if r_par.dangerous_prc
            w_base.w_prc = W_PRC();
            w_base.w_prc.alpha_star = 0;
            w_base.w_prc.vus = 0.5;
        else
            compute_prc_quantities( h_c, a_par, w_base );
        end
    end
    
    %
    % Plotting
    %
    
    if a_par.draw_scatter_plots
        draw_scatter_plot( a_par, r_par, h_c);
        if (w_base.w_prc ~= 0)
            draw_prc_plot(a_par, r_par, w_base.w_prc);
        end
        if (w_base.w_bc ~= 0)
            draw_roc_plot(a_par, r_par, w_base.w_bc);
        end
    end
    
    %
    % cleanup
    %
    
    %work_list{r_par.run_number} = w_base;
    
    
    
    switch a_par.evaluation_mode
        case 'train'
            switch a_par.opt_value
                case 'alpha_star', val = -w_base.w_prc.alpha_star;
                case 'vus',    val = -w_base.w_prc.vus;
                case 'ter',    val = -w_base.w_bc.best_ter;
                case 'ber',    val = -w_base.w_bc.best_ber;
                case 'f_score',val = -w_base.w_bc.best_f_score;

                otherwise
                    warning('%s not recognized!', a_par.opt_value)
            end
            
        case 'test'
            val = zeros(5, 1);
            val(1) = w_base.w_prc.alpha_star;
            val(2) = w_base.w_prc.vus;
            val(3) = w_base.w_bc.best_ter;
            val(4) = w_base.w_bc.best_ber;
            val(5) = w_base.w_bc.best_f_score;
            
        otherwise
            warning('%s not recognized!\n', a_par.evaluation_mode)
    end
    
    if a_par.save_workspace
        save('w_base.mat', 'w_base');
    end
    
    
    fprintf('Run %s done after %0.2f seconds.\n', r_par.title_tag, toc);


end

