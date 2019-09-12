classdef Histogram_Container < handle
    %HISTOGRAM_CONTAINER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        P_hist_A;
        P_hist_B;
        P_hist_N;
        
        R_hist_A;
        R_hist_B;
        R_hist_N;
        
        F_hist_A;
        F_hist_N;
        
        A_thresholds;
        B_thresholds;
    end
    
    methods
        function [ h_c ] = Histogram_Container(a_par, r_par)
            h_c.P_hist_A = [];
            h_c.P_hist_B = [];
            h_c.P_hist_N = [];
            
            h_c.R_hist_A = [];
            h_c.R_hist_B = [];
            h_c.R_hist_N = [];
            
            h_c.F_hist_A = [];
            h_c.F_hist_N = [];
            
            tic;
            fprintf('Beginning histogram construction.\n');
            
            switch a_par.model
                case 'mmt'
                    coeff = r_par.pred_fun.gabor_coeff;
                case 'kolm'
                    coeff = r_par.coeffs;
            end
            
            
            r = (size(coeff, 1)-1)/2;
    
            %L = r_par.gabor_L_steps;
            %T = a_par.gap_T_steps;

            switch a_par.evaluation_mode 
                case 'train'
                    filenums = a_par.training_file_set;
                case 'test'
                    filenums = a_par.testing_file_set;
                otherwise
                    warning('%s not reocgnized!\n', a_par.evaluation_mode)
            end
            
            switch a_par.model
                case 'mmt'
                    for m = 1:length(filenums)
                
                        cur_filenum = filenums(m);
                        filename = sprintf('%s/sim_%d.dat', ...
                            a_par.datapath, cur_filenum);
                        fprintf('Loading file %s (%d out of %d).\n', ...
                            filename, m, length(filenums));
                        %YY;
                        load(filename, 'YY', '-mat');

                        B = nan(size(YY));

                        TT = h_c.get_T_eval_indices(a_par);

                        for t = TT
                            cur_slice = YY(:, t);
                            repslice = [cur_slice((end-r+1):end); cur_slice; ...
                                        cur_slice(1:r)];
                            B(:, t) = conv(repslice, coeff, 'valid');
                        end

                        B = abs(B);
                        A = abs(YY);
                        clear YY

                        % Doesn't do anything!
                        %MM = h_c.construct_mask(A, B, a_par, r_par);

                        %
                        % 
                        %

                        build_hist_max_exact(h_c, A, B, a_par, r_par);
                        clear A B YY;

                    end
                    
                case 'kolm'
                    for m = 1:length(filenums)
                
                        cur_filenum = filenums(m);
                        filename = sprintf('%s/run_%d.dat', ...
                            a_par.datapath, cur_filenum);
                        fprintf('Loading file %s (%d out of %d).\n', ...
                            filename, m, length(filenums));
                        load(filename, 'e_diss_mu', 'f_u12', '-mat');
                        
                        B = zeros(size(e_diss_mu));
                        
                        for t = 1:length(e_diss_mu)
                            switch a_par.fourier_set
                                case 'full'
                                    warning('This doesn''t work!\n')
                                case 'half'
                                    cur_f = f_u12(t, :, :);
                                case 'quadrant'
                                    L = size(f_u12, 2);
                                    LL =  ceil(L/2):L;
                                    cur_f = f_u12(t, LL, :);
                            end
                            
                            B(t) = sum(cur_f(:).*coeff');
                        end
                        
                        B = abs(B);
                        A = abs(e_diss_mu);
                        clear e_diss_mu f_u12
                        
                        build_hist_max_exact(h_c, A, B, a_par, r_par);
                        clear A B YY;
                    end
            end
            
            
            
            switch a_par.threshold_distribution
                case 'linear'
                    min_a = min(h_c.P_hist_A(:));
                    max_a = max(h_c.P_hist_A(:));
                    min_b = min(h_c.P_hist_B(:));
                    max_b = max(h_c.P_hist_B(:));

                    h_c.A_thresholds = linspace(min_a, max_a, a_par.n_a_thresholds);
                    h_c.B_thresholds = linspace(min_b, max_b, a_par.n_b_thresholds);

                otherwise
                    warning('%s not recognized!', a_par.threshold_distribution)
            end
            
            fprintf('Histogram construction done after %0.2f seconds.\n', toc)
            
            
        end
        
        
        
        function [ outcode ] = build_hist_max_exact(h_c, A, B, a_par, r_par)
            %size_unrolled = (6000 - a_par.drop_T_steps + 1);
            
            switch a_par.model
                case 'mmt'
                    max_A = max(A, [], 1);
                    max_B = max(B, [], 1);
                    
                    % for which spatial max is superfluous anyway...
                case 'kolm'
                    max_A = A;
                    max_B = B;
            end
            
            %t_start = a_par.drop_T_steps + 1;
            %t_end = 6000;
            TT = h_c.get_T_eval_indices(a_par);
            

            size_unrolled = (a_par.data_length - a_par.drop_T_steps + 1);
  

            corr_B = nan(size_unrolled, 1);
            corr_A = nan(size_unrolled, 1);

            
            %for t = t_start:(t_end - a_par.gap_T_steps)                
            %for t = TT(1:(end-a_par.gap_T_steps))
            for index = 1:(length(TT)-a_par.gap_T_steps)
                %index = (t - t_start) + 1;
                
                corr_A(index) = max_A(TT(index) + a_par.gap_T_steps);
                corr_B(index) = max_B(TT(index));
            end
            
%             figure(22);
%             clf;
%             hold on
%             plot(max_A);
%             plot(max_B);
%             hold off
%             
%             figure(23);
%             clf;
%             histogram2(corr_A, corr_B);
            
            % P and R histogram
            
            %[N, A_edges, B_edges] = histcounts2(corr_A,corr_B,a_par.data_hist_bins, 'BinMethod', 'auto');
            [N, A_edges, B_edges] = histcounts2(corr_A,corr_B,a_par.data_hist_bins);
            
            AA  = zeros(length(N(:)), 1);
            BB  = zeros(length(N(:)), 1);
            NN  = zeros(length(N(:)), 1);

            AA_center = (A_edges(2:end) + A_edges(1:end-1))/2;
            BB_center = (B_edges(2:end) + B_edges(1:end-1))/2;

            for k = 1:size(N, 1)
                for j = 1:size(N, 2)
                    q = k + (j-1)*size(N, 1);
                    AA(q) = AA_center(k);
                    BB(q) = BB_center(j);
                    NN(q) = N(k,j);
                end
            end

            is_filled = (NN > 0);
            AA = AA(is_filled);
            BB = BB(is_filled);
            NN = NN(is_filled);
            
            
            h_c.P_hist_A = [h_c.P_hist_A; AA];
            h_c.P_hist_B = [h_c.P_hist_B; BB];
            h_c.P_hist_N = [h_c.P_hist_N; NN];
            
            %h_c.R_hist_A = [h_c.R_hist_A; AA];
            %h_c.R_hist_B = [h_c.R_hist_B; BB];
            %h_c.R_hist_N = [h_c.R_hist_N; NN];
            
            clear AA BB NN is_filled;
            
%             % F histogram
%             
%             %[N, A_edges] = histcounts(corr_A,a_par.data_hist_bins, 'BinMethod', 'auto');
%             [N, A_edges] = histcounts(corr_A,a_par.data_hist_bins);
%             
%             AA  = zeros(length(N(:)), 1);
%             NN  = zeros(length(N(:)), 1);
%             
%             AA_center = (A_edges(2:end) + A_edges(1:end-1))/2;
%             
%             for k = 1:length(N)
%                 AA(k) = AA_center(k);
%                 NN(k) = N(k);
%             end
% 
%             is_filled = (NN > 0);
%             AA = AA(is_filled);
%             NN = NN(is_filled);
%             
%             h_c.F_hist_A = [h_c.F_hist_A; AA];
%             h_c.F_hist_N = [h_c.F_hist_N; NN];
            
            clear corr_A corr_B
            clear max_A max_B
            
            outcode = 1;
        end
        
        
        function [ tp, tn, fp, fn ] = get_binary_classification(h_c, a_par, a, b, hist_type)
            if a_par.use_identical_histograms
                hist_type = 'p';
            end
            
            switch hist_type
                case 'p'
                    TP_mask = (h_c.P_hist_B >= b) & (h_c.P_hist_A >= a);
                    TN_mask = (h_c.P_hist_B < b) & (h_c.P_hist_A < a);
                    FP_mask = (h_c.P_hist_B >= b) & (h_c.P_hist_A < a);
                    FN_mask = (h_c.P_hist_B < b) & (h_c.P_hist_A >= a);
                    tp = sum(h_c.P_hist_N(TP_mask));
                    tn = sum(h_c.P_hist_N(TN_mask));
                    fp = sum(h_c.P_hist_N(FP_mask));
                    fn = sum(h_c.P_hist_N(FN_mask));
                    
                case 'r'
                    TP_mask = (h_c.R_hist_B >= b) & (h_c.R_hist_A >= a);
                    TN_mask = (h_c.R_hist_B < b) & (h_c.R_hist_A < a);
                    FP_mask = (h_c.R_hist_B >= b) & (h_c.R_hist_A < a);
                    FN_mask = (h_c.R_hist_B < b) & (h_c.R_hist_A >= a);
                    tp = sum(h_c.R_hist_N(TP_mask));
                    tn = sum(h_c.R_hist_N(TN_mask));
                    fp = sum(h_c.R_hist_N(FP_mask));
                    fn = sum(h_c.R_hist_N(FN_mask));
                    
                case 'f'
                    TP_mask = (h_c.F_hist_A >= a);
                    TN_mask = (h_c.F_hist_A < a);
                    tp = sum(h_c.F_hist_N(TP_mask));
                    tn = sum(h_c.F_hist_N(TN_mask));
                    
                otherwise
                    warning('%s not recognized!\n', hist_type)
            end
        end
        
        
        
        function [ vals ] = get_binary_classification_full_set(h_c, a_par)
            if ~a_par.use_identical_histograms
                warning('This isn''t implemented!\n');
            end
            
            vals = struct;
            vals.a_thresh = h_c.A_thresholds;
            vals.b_thresh = h_c.B_thresholds;
            vals.tp = zeros(length(h_c.A_thresholds), length(h_c.B_thresholds));
            vals.tn = zeros(length(h_c.A_thresholds), length(h_c.B_thresholds));
            vals.fp = zeros(length(h_c.A_thresholds), length(h_c.B_thresholds));
            vals.fn = zeros(length(h_c.A_thresholds), length(h_c.B_thresholds));
            
            for k_data = 1:length(h_c.P_hist_A)
                cur_a = h_c.P_hist_A(k_data);
                cur_b = h_c.P_hist_B(k_data);
                cur_n = h_c.P_hist_N(k_data);
                
                TP_mask = (cur_b >= h_c.B_thresholds) & (cur_a >= h_c.A_thresholds)';
                TN_mask = (cur_b < h_c.B_thresholds) & (cur_a < h_c.A_thresholds)';
                FP_mask = (cur_b >= h_c.B_thresholds) & (cur_a < h_c.A_thresholds)';
                FN_mask = (cur_b < h_c.B_thresholds) & (cur_a >= h_c.A_thresholds)';
                
                vals.tp(TP_mask) = vals.tp(TP_mask) + cur_n;
                vals.tn(TN_mask) = vals.tn(TN_mask) + cur_n;
                vals.fp(FP_mask) = vals.fp(FP_mask) + cur_n;
                vals.fn(FN_mask) = vals.fn(FN_mask) + cur_n;
            end

            
        end
        
        
        function [TT] = get_T_eval_indices(h_c, a_par)
            TT = [];
            switch a_par.domain_frac
                case 'full'
                    t_start = a_par.drop_T_steps + 1;
                    t_end = a_par.data_length;
                    TT = t_start:t_end;
                case 'half'
                    t_start = a_par.drop_T_steps + 1;
                    t_end = floor(a_par.data_length/2) + floor(t_start/2);
                    TT = t_start:t_end;
                otherwise
                    warning('%s not recognized!\n', a_par.domain_frac);
            end

        end
    end
    
end

