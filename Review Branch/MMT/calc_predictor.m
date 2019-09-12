function [B] = calc_predictor(ab_par, YY)
%CALC_AB Summary of this function goes here
%   Detailed explanation goes here

    B = nan(size(YY));
    r = (size(ab_par.gabor_coefficients, 1)-1)/2;
    
    for t = ab_par.t_start:ab_par.t_end
        cur_slice = YY(:, t);
        repslice = [cur_slice((end-r+1):end); cur_slice; cur_slice(1:r)];
        B(:, t) = conv(repslice, ab_par.gabor_coefficients, 'valid');
    end

    B = abs(B);
end

