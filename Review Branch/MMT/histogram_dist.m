function [ AA,  NN] = histogram_dist( h_par, A)
%HISTOGRAM_IND_PRE Summary of this function goes here
%   Detailed explanation goes here

    
    [N, A_edges] = histcounts(A,h_par.n_bins, 'BinMethod', 'auto');
    
    
    
    AA  = zeros(length(N(:)), 1);
    NN  = zeros(length(N(:)), 1);
    


    AA_center = (A_edges(2:end) + A_edges(1:end-1))/2;
    
    
    for k = 1:length(N)
        AA(k) = AA_center(k);
        NN(k) = N(k);
    end
    
    is_filled = (NN > 0);
    AA = AA(is_filled);
    NN = NN(is_filled);

end

