function [outcode] = compute_bc_quantities( h_c, a_par, w_base )
%COMPUTE_BC_QUANTITIES Summary of this function goes here
%   Detailed explanation goes here

    w_bc = W_BC();

    w_bc.A_threshold = a_par.fixed_A_threshold;
    w_bc.B_threshold = h_c.B_thresholds;
    a_hat = a_par.fixed_A_threshold;
    
    for k = 1:a_par.n_b_thresholds
        b_hat = w_bc.B_threshold(k);
        [ tp, tn, fp, fn ] = h_c.get_binary_classification(a_par, a_hat, b_hat, 'p');
        w_bc.p_true_positive(k) = tp;
        w_bc.p_false_positive(k)= fp;
        w_bc.p_true_negative(k) = tn;
        w_bc.p_false_negative(k)= fn;
        
        [ tp, tn, fp, fn ] = h_c.get_binary_classification(a_par, a_hat, b_hat, 'r');
        w_bc.r_true_positive(k) = tp;
        w_bc.r_false_positive(k)= fp;
        w_bc.r_true_negative(k) = tn;
        w_bc.r_false_negative(k)= fn;
    end
    
    w_bc.precision = w_bc.p_true_positive ./ ...
        (w_bc.p_true_positive +  w_bc.p_false_positive);
    w_bc.recall = w_bc.r_true_positive ./ ...
        (w_bc.r_true_positive +  w_bc.r_false_negative);
    w_bc.specificity = w_bc.p_true_negative ./ ...
        (w_bc.p_true_negative +  w_bc.p_false_positive);
    w_bc.total_error_rate = (w_bc.p_true_positive + w_bc.p_true_negative)./ ...
        (w_bc.p_true_negative +  w_bc.p_false_positive + ...
        w_bc.r_true_positive +  w_bc.r_false_negative);
    
    w_bc.precision(isnan(w_bc.precision)) = 0;
    w_bc.recall(isnan(w_bc.recall)) = 0;    
    w_bc.specificity(isnan(w_bc.specificity)) = 0;
    w_bc.total_error_rate(isnan(w_bc.total_error_rate)) = 0;
    
    w_bc.balanced_error_rate = (w_bc.specificity+w_bc.recall)/2;
    w_bc.F_score = (1 + a_par.F_score_beta^2)*(w_bc.precision.*w_bc.recall)./ ...
        (a_par.F_score_beta^2 * w_bc.precision + w_bc.recall);
    
    w_bc.best_ter = max(w_bc.total_error_rate(:));
    w_bc.best_ber = max(w_bc.balanced_error_rate(:));
    w_bc.best_f_score = max(w_bc.F_score(:));
    
    w_base.w_bc = w_bc;
    
    if (a_par.verbosity >= 1)
        fprintf('Final TER = %0.3f.\n', ...
            w_bc.best_ter);
        fprintf('Final BER = %0.3f.\n', ...
            w_bc.best_ber);
        fprintf('Final F_beta = %0.3f (beta = %0.1f).\n', ...
            w_bc.best_f_score, a_par.F_score_beta);
        fprintf('BC calculations done after %0.2f seconds.\n', toc);
    end
    
    
    outcode = 1;
end

