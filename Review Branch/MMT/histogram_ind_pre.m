function [ AA, BB, NN] = histogram_ind_pre( h_par, A, B)
%HISTOGRAM_IND_PRE Summary of this function goes here
%   Detailed explanation goes here

    
    [N, A_edges, B_edges] = histcounts2(A,B,h_par.n_bins, 'BinMethod', 'auto');
    
    
    
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

end

