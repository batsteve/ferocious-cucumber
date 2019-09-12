function [ outcode ] = compute_basic_statistics( h_c, a_par, w_base )
%COMPUTE_BASIC_STATISTICS Summary of this function goes here
%   Detailed explanation goes here

    %if (nargin == 4)
    %    NN = ones(size(AA));
    %end
    
    if (a_par.verbosity >= 1)
        tic;
        fprintf('Starting basic statistics calculations.\n');
    end
    
    AA = h_c.P_hist_A;
    BB = h_c.P_hist_B;
    NN = h_c.P_hist_N;

    w_stat = W_Statistics();
    w_base.w_stat = w_stat;
    
    w_stat.full_counts = sum(NN);
    w_stat.full_B_mean = sum(BB.*NN)/w_stat.full_counts;
    w_stat.full_B_var = sum(((BB - w_stat.full_B_mean).^2.*NN))/w_stat.full_counts;
    w_stat.full_A_mean = sum(AA.*NN)/w_stat.full_counts;
    w_stat.full_A_var = sum(((AA - w_stat.full_A_mean).^2.*NN))/w_stat.full_counts;
    
    if (a_par.verbosity >= 1)
        fprintf('Basic statistics calculations done after %0.2f seconds.\n', toc);
    end
    
    outcode = 1;

end

