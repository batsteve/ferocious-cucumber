function [ AA, LL ] = parametersToWavelets( a_par, coeffs)
%PARAMETERSTOWAVELETS Summary of this function goes here
%   Detailed explanation goes here

    num_wavelets = a_par.hypothesis_class_size;
    
    AA = zeros(num_wavelets, 1);
    AA(2:end) = coeffs(1:(num_wavelets-1));
    AA(1) = real(sqrt(1 - sum(AA(:)).^2));
    
    CC = zeros(num_wavelets, 1);
    
    switch a_par.parameter_unroll_coding
        case 'independent'
            CC =  coeffs(num_wavelets:end);
        case 'stacked'
            for k = 1:num_wavelets
                CC(k) = sum(coeffs(num_wavelets:(num_wavelets + k - 1)));
            end
            
            k_off = log(0.5);
            B = 0.2;
            A_scale = log(2000) - log(0.5);
            
            CC = (CC)/(CC(num_wavelets) + B)*A_scale + k_off;
            
        otherwise
            warning('%s unrecognized!\n', a_par.parameter_unroll_coding)
    end
    
    
    
    switch a_par.parameter_scale
        case 'log'
            LL = exp(CC);
        case 'linear'
            LL = CC;
        otherwise
            warning('%s not recognized!\n', a_par.parameter_scale);
    end
end

