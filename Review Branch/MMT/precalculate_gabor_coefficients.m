function [ outcode ] = precalculate_gabor_coefficients( ab_par )
%PRECALCULATE_GBUR_COEFFICIENTS Summary of this function goes here
%   Detailed explanation goes here

    %x_d = min(20, ceil(par.gbur_L*4));     %This was embarassing.
    x_d = min(floor(ab_par.n_x/2), max(20, ceil(ab_par.gabor_L*4)));
    t_d = 0;
    
    xx = linspace(-x_d, x_d, 2*x_d +1);
    tt = linspace(-t_d, t_d, 2*t_d +1);
    
    xx = min(abs(xx), ab_par.n_x - abs(xx));
    
    [gbur_corr_x, gbur_corr_t] = ndgrid(xx, tt);
    
    switch ab_par.predictor_kernel
        case 'sinc'
            G = @(x, t) sinc(x/ab_par.gabor_L);
        case 'exp_squared'
            G = @(x, t) exp(-x.^2./(2*ab_par.gabor_L^2));
    end
    
    % There was a problem with imaginary numbers, somewhere?
    ab_par.gabor_coefficients = G(gbur_corr_x, gbur_corr_t);
    
    gabor_norm = sum(ab_par.gabor_coefficients(:).^2);
    ab_par.gabor_coefficients = ab_par.gabor_coefficients./gabor_norm;
    
    outcode = 1;

end

