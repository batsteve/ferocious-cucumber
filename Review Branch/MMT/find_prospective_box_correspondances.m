function [ A_from_B, B_from_B ] = find_prospective_box_correspondances(b_par, A, B)
%FIND_BOX_CORRESPONDANCES Summary of this function goes here
%   Detailed explanation goes here

    size_unrolled = (b_par.t_end - b_par.t_start + 1)*b_par.n_x;

    B_from_B = nan(size_unrolled, 1);
    A_from_B = nan(size_unrolled, 1);
    
    RR = -b_par.big_delta_x:b_par.big_delta_x;
    XX = 1:b_par.n_x;

    for t = b_par.t_start:(b_par.t_end - b_par.big_delta_t_min)
        A_slice = A(:, t + b_par.big_delta_t_min);
        A_slice = [A_slice, A_slice, A_slice];
        
        index = XX + (t - b_par.t_start)*b_par.n_x;
        A_from_B(index) = max(A_slice(b_par.n_x + XX + RR'));
        B_from_B(index) = B(XX, t);
    end


end

