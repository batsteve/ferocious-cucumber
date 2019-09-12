function [ A_from_A, B_from_A ] = find_retrospective_box_correspondances(b_par, A, B)
%FIND_BOX_CORRESPONDANCES Summary of this function goes here
%   Detailed explanation goes here

    size_unrolled = (b_par.t_end - b_par.t_start + 1)*b_par.n_x;

    B_from_A = nan(size_unrolled, 1);
    A_from_A = nan(size_unrolled, 1);
    
    RR = -b_par.big_delta_x:b_par.big_delta_x;
    XX = 1:b_par.n_x;

    for t = (b_par.t_start + b_par.big_delta_t_min):b_par.t_end
        
        B_slice = B(:, t - b_par.big_delta_t_min);
        B_slice = [B_slice, B_slice, B_slice];
        
        
        index = XX + (t - (b_par.t_start + b_par.big_delta_t_min))*b_par.n_x;
        A_from_A(index) = A(XX, t);
        B_from_A(index) = max(B_slice(b_par.n_x + XX + RR'));   
    end




end

