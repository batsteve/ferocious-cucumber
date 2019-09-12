function [ outcode ] = compute_prc_quantities( h_c, a_par, w_base )
%COMPUTE_PRC_QUANTITIES Summary of this function goes here
%   Detailed explanation goes here

    %
    % Declaration of constants
    %

    %guess_auc = 0.5;

    if (a_par.verbosity >= 1)
        tic;
        fprintf('Starting PRC calculations.\n');
    end
    
    %total_events = sum(h_c.F_hist_N(:));

    
    w_prc = W_PRC();
    w_prc.predictor_polarity = 'positive';
    

    %
    % check every pair of thresholds
    %
    
    switch a_par.bin_count_method
        case 'one-by-one'
            
            w_prc.A_threshold = h_c.A_thresholds;
            w_prc.B_threshold = h_c.B_thresholds;



            w_prc.p_true_positive = zeros(a_par.n_a_thresholds, a_par.n_b_thresholds);
            w_prc.p_true_negative = zeros(a_par.n_a_thresholds, a_par.n_b_thresholds);
            w_prc.p_false_positive = zeros(a_par.n_a_thresholds, a_par.n_b_thresholds);
            w_prc.p_false_negative = zeros(a_par.n_a_thresholds, a_par.n_b_thresholds);

            w_prc.r_true_positive = zeros(a_par.n_a_thresholds, a_par.n_b_thresholds);
            w_prc.r_true_negative = zeros(a_par.n_a_thresholds, a_par.n_b_thresholds);
            w_prc.r_false_positive = zeros(a_par.n_a_thresholds, a_par.n_b_thresholds);
            w_prc.r_false_negative = zeros(a_par.n_a_thresholds, a_par.n_b_thresholds);

            w_prc.f_positive = zeros(a_par.n_a_thresholds, a_par.n_b_thresholds);
            w_prc.f_negative = zeros(a_par.n_a_thresholds, a_par.n_b_thresholds);
    
            for k = 1:a_par.n_a_thresholds
                for j = 1:a_par.n_b_thresholds

                    a_hat = w_prc.A_threshold(k);
                    b_hat = w_prc.B_threshold(j);
                    [ tp, tn, fp, fn ] = h_c.get_binary_classification(a_par, a_hat, b_hat, 'p');
                    w_prc.p_true_positive(k, j) = tp;
                    w_prc.p_false_positive(k, j)= fp;
                    w_prc.p_true_negative(k, j) = tn;
                    w_prc.p_false_negative(k, j)= fn;

                    [ tp, tn, fp, fn ] = h_c.get_binary_classification(a_par, a_hat, b_hat, 'r');
                    w_prc.r_true_positive(k, j) = tp;
                    w_prc.r_false_positive(k, j)= fp;
                    w_prc.r_true_negative(k, j) = tn;
                    w_prc.r_false_negative(k, j)= fn;

                    [ tp, tn, fp, fn ] = h_c.get_binary_classification(a_par, a_hat, b_hat, 'f');
                    w_prc.f_positive(k, j) = (tp + fn);
                    w_prc.f_negative(k, j) = (tn + fp);
                end
            end
    
        case 'full_set'
            vals = h_c.get_binary_classification_full_set(a_par);
            
            w_prc.A_threshold = vals.a_thresh;
            w_prc.B_threshold = vals.b_thresh;
            
            w_prc.p_true_positive = vals.tp;
            w_prc.p_false_positive= vals.fp;
            w_prc.p_true_negative = vals.tn;
            w_prc.p_false_negative= vals.fn;
            
            w_prc.r_true_positive = vals.tp;
            w_prc.r_false_positive= vals.fp;
            w_prc.r_true_negative = vals.tn;
            w_prc.r_false_negative= vals.fn;
            
            w_prc.f_positive = (vals.tp + vals.fn);
            w_prc.f_negative = (vals.tn + vals.fp);
        otherwise
            warning('Not recognized!\n')
    end
    
    
    
    w_prc.precision = zeros(a_par.n_a_thresholds, a_par.n_b_thresholds);
    w_prc.recall = zeros(a_par.n_a_thresholds, a_par.n_b_thresholds);
    w_prc.ferocity = zeros(a_par.n_a_thresholds, a_par.n_b_thresholds);
    w_prc.specificity = zeros(a_par.n_a_thresholds, a_par.n_b_thresholds);
    w_prc.pre_ferocity = zeros(a_par.n_a_thresholds, a_par.n_b_thresholds);
    w_prc.dangerosity = zeros(a_par.n_a_thresholds, a_par.n_b_thresholds);



    %
    % Calc derived quantities
    %
    
    w_prc.precision = w_prc.p_true_positive ./ ...
        (w_prc.p_true_positive +  w_prc.p_false_positive);
    w_prc.recall = w_prc.r_true_positive ./ ...
        (w_prc.r_true_positive +  w_prc.r_false_negative);
    w_prc.ferocity = (w_prc.f_positive) ./ ...
        (w_prc.f_positive +  w_prc.f_negative);
    
    
    w_prc.specificity = w_prc.p_true_negative ./ ...
        (w_prc.p_true_negative +  w_prc.p_false_positive);
    w_prc.pre_ferocity = (w_prc.p_true_positive + w_prc.p_false_positive) ./ ...
        (w_prc.p_true_positive +  w_prc.p_false_positive + ...
        w_prc.p_true_negative +  w_prc.p_false_negative);
    w_prc.dangerosity = (w_prc.p_true_positive + w_prc.p_false_negative) ./ ...
        (w_prc.p_true_positive +  w_prc.p_false_positive +  w_prc.p_false_negative);
    
    w_prc.precision(isnan(w_prc.precision)) = 0;
    w_prc.recall(isnan(w_prc.recall)) = 0;
    w_prc.ferocity(isnan(w_prc.ferocity)) = 0;
    
    w_prc.specificity(isnan(w_prc.specificity)) = 0;
    w_prc.pre_ferocity(isnan(w_prc.pre_ferocity)) = 0;
    w_prc.dangerosity(isnan(w_prc.dangerosity)) = 0;
    
    if (a_par.verbosity >= 2)
        fprintf('PRC threshold pair sums done after %0.2f seconds.\n', toc);
    end
    
    
    
    switch a_par.auc_integration_algorithm
        case 'scattered interpolant'
            %
            % Ferocity auc
            %
            
            R = linspace(0, 1, a_par.auc_integration_points);
            F = linspace(0, 1, a_par.auc_integration_points);
            
            [RR, FF] = meshgrid(R, F);
            data_R = w_prc.recall(:);
            data_F = w_prc.ferocity(:);
            data_P = w_prc.precision(:);
            II = scatteredInterpolant(data_R, data_F, data_P);
            PP = II(RR, FF);
            
            PP(PP > 1) = 1;
            PP(PP < 0) = 0;
            
            % cuz edge points and stuff
            
            center_PP = PP(2:(end-1), 2:(end-1));
            side_PP = [PP(2:(end-1), 1); PP(2:(end-1), end); ...
                PP(1, 2:(end-1))';  PP(end, 2:(end-1))'];
            corner_PP = [PP(1, 1); PP(1, end); PP(end, 1);  PP(end, end)];
            
            w_prc.vus = (sum(center_PP(:)) + 1/2*sum(side_PP(:)) + 1/4*sum(corner_PP(:))) /...
                a_par.auc_integration_points^2;
              
            center_PP = PP(:, 2:(end-1));
            side_PP = [PP(:, 1), PP(:, end)];
            
            alpha_list = zeros(a_par.auc_integration_points, 1);
            search_KK = 2:(a_par.auc_integration_points);
            % don't check fer=0, because of 0/0 numerical stuff
            for k = search_KK
                alpha_list(k) = (sum(center_PP(k, :), 2) + 1/2*sum(side_PP(k, :), 2))./a_par.auc_integration_points - F(k);
            end
            
            [w_prc.alpha_star, pre] = max(alpha_list(:));
            fprintf('a_star = %0.3f, at index %d, fer=%0.3f.\n', ...
                w_prc.alpha_star, pre, F(pre));
            
              
            
        case 'trapezoid'
            %weights = [w_prc.recall(:, 1)/2, ...
            %    w_prc.recall(:, 3:end) - w_prc.recall(:, 1:end-2),...
            %    w_prc.recall(:, end)/2]';
            diffs = diff(w_prc.recall, 1, 2);
            weights = -([diffs, zeros(size(diffs, 1), 1)] + ...
                 [zeros(size(diffs, 1), 1), diffs])/2;
            F_alpha = sum(weights.*w_prc.precision, 2);
            
            v_diffs = diff(w_prc.ferocity, 1, 1);
            v_weights = -([v_diffs; zeros(1, size(v_diffs, 2))] + ...
                 [zeros(1, size(v_diffs, 2)); v_diffs])/2;
            w_prc.vus = sum(v_weights(1, :)'.*F_alpha);
            
            alpha_list = zeros(a_par.n_a_thresholds, 1);
            search_KK = 2:(a_par.n_a_thresholds - 1);
            for k = search_KK
                alpha_list(k) = F_alpha(k) - w_prc.ferocity(k, 1);
            end
            
            [w_prc.alpha_star, pre] = max(alpha_list(:));
            fprintf('a_star = %0.3f, at index %d, fer=%0.3f.\n', ...
                w_prc.alpha_star, pre, w_prc.ferocity(pre, 1));
  
        otherwise
            warning('%s not recognized!', a_par.auc_integration_algorithm);
    end
    
    
    
    if (w_prc.alpha_star > 1)
        warning('Alpha* is too high!\n');
        disp(PP)
    end
    
            
    
    w_base.w_prc = w_prc;
    
    if (a_par.verbosity >= 1)
        fprintf('Final VUS = %0.3f (%s polarity).\n', ...
            w_prc.vus, w_prc.predictor_polarity);
        fprintf('Final alpha* = %0.3f (%s polarity).\n', ...
            w_prc.alpha_star, w_prc.predictor_polarity);
        fprintf('PRC calculations done after %0.2f seconds.\n', toc);
    end
    
    outcode = 1;

end







            
%             %
%             % count the events
%             % mask rule depends on predictor polarity
%             %
%             
%             switch w_prc.predictor_polarity
%                 case 'positive'
%                     TP_mask = (h_c.P_hist_B >= w_prc.B_threshold(j)) & (h_c.P_hist_A >= w_prc.A_threshold(k));
%                     TN_mask = (h_c.P_hist_B < w_prc.B_threshold(j)) & (h_c.P_hist_A < w_prc.A_threshold(k));
%                     FP_mask = (h_c.P_hist_B >= w_prc.B_threshold(j)) & (h_c.P_hist_A < w_prc.A_threshold(k));
%                     FN_mask = (h_c.P_hist_B < w_prc.B_threshold(j)) & (h_c.P_hist_A >= w_prc.A_threshold(k));
%                 case 'negative'
%                     TP_mask = (h_c.P_hist_B < w_prc.B_threshold(j)) & (h_c.P_hist_A >= w_prc.A_threshold(k));
%                     TN_mask = (h_c.P_hist_B >= w_prc.B_threshold(j)) & (h_c.P_hist_A < w_prc.A_threshold(k));
%                     FP_mask = (h_c.P_hist_B < w_prc.B_threshold(j)) & (h_c.P_hist_A < w_prc.A_threshold(k));
%                     FN_mask = (h_c.P_hist_B >= w_prc.B_threshold(j)) & (h_c.P_hist_A >= w_prc.A_threshold(k));
%                 otherwise
%                     warning('%s not recognized!', w_prc.predictor_polarity);
%             end
%             
%             w_prc.p_true_positive(k, j) = sum(h_c.P_hist_N(TP_mask));
%             w_prc.p_false_positive(k, j)= sum(h_c.P_hist_N(FP_mask));
%             w_prc.p_true_negative(k, j) = sum(h_c.P_hist_N(TN_mask));
%             w_prc.p_false_negative(k, j)= sum(h_c.P_hist_N(FN_mask));
%             
%             
%             
%             switch w_prc.predictor_polarity
%                 case 'positive'
%                     TP_mask = (h_c.R_hist_B >= w_prc.B_threshold(j)) & (h_c.R_hist_A >= w_prc.A_threshold(k));
%                     TN_mask = (h_c.R_hist_B < w_prc.B_threshold(j)) & (h_c.R_hist_A < w_prc.A_threshold(k));
%                     FP_mask = (h_c.R_hist_B >= w_prc.B_threshold(j)) & (h_c.R_hist_A < w_prc.A_threshold(k));
%                     FN_mask = (h_c.R_hist_B < w_prc.B_threshold(j)) & (h_c.R_hist_A >= w_prc.A_threshold(k));
%                 case 'negative'
%                     TP_mask = (h_c.R_hist_B < w_prc.B_threshold(j)) & (h_c.R_hist_A >= w_prc.A_threshold(k));
%                     TN_mask = (h_c.R_hist_B >= w_prc.B_threshold(j)) & (h_c.R_hist_A < w_prc.A_threshold(k));
%                     FP_mask = (h_c.R_hist_B < w_prc.B_threshold(j)) & (h_c.R_hist_A < w_prc.A_threshold(k));
%                     FN_mask = (h_c.R_hist_B >= w_prc.B_threshold(j)) & (h_c.R_hist_A >= w_prc.A_threshold(k));
%                 otherwise
%                     warning('%s not recognized!', w_prc.predictor_polarity);
%             end
%             
%             w_prc.r_true_positive(k, j) = sum(h_c.R_hist_N(TP_mask));
%             w_prc.r_false_positive(k, j)= sum(h_c.R_hist_N(FP_mask));
%             w_prc.r_true_negative(k, j) = sum(h_c.R_hist_N(TN_mask));
%             w_prc.r_false_negative(k, j)= sum(h_c.R_hist_N(FN_mask));
%             
%             P_mask = (h_c.F_hist_A >= w_prc.A_threshold(k));
%             N_mask = (h_c.F_hist_A < w_prc.A_threshold(k));
% 
%             w_prc.f_positive(k, j) = sum(h_c.F_hist_N(P_mask));
%             w_prc.f_negative(k, j) = sum(h_c.F_hist_N(N_mask));