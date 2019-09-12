classdef Predictor_Container < handle
    %PREDICTOR_CONTAINER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        coefficient_list;
        gabor_coeff;
        
        % coefficient_list[1] = alpha
        % coefficient_list[2] = L_1
        % coefficient_list[3] = L_2
        
        % G = (1-alpha) G[L_1] + alpha G[L_2]
    end
    
    methods
        function [ p_c ] = Predictor_Container(a_par, new_coeffs)
            p_c.coefficient_list = new_coeffs;
            
            p_c.gabor_coeff = p_c.precalculate_gabor_coefficients(a_par);
        end
        
        function [ gabor_coefficients ] = precalculate_gabor_coefficients( p_c, a_par )
        %PRECALCULATE_GABOR_COEFFICIENTS Summary of this function goes here
        %   Detailed explanation goes here
        
            num_wavelets = a_par.hypothesis_class_size;
            
            [ AA, LL ] = parametersToWavelets( a_par, p_c.coefficient_list );
            
            max_L = max(LL(:));

            x_d = min(floor(a_par.n_x/2), max(20, ceil(max_L*4)));
            xx = linspace(-x_d, x_d, 2*x_d +1);
            xx = min(abs(xx), a_par.n_x - abs(xx));

            wavelets = zeros(length(xx), num_wavelets);
            for k = 1:num_wavelets

                switch a_par.kernel_shape
                    case 'sinc'
                        G = @(x) sinc(x/LL(k));
                    case 'exp_squared'
                        G = @(x) exp(-x.^2./(2*LL(k)^2));
                    otherwise
                        warning('%s not recognized!');
                end

                % There was a problem with imaginary numbers, somewhere?
                raw_wavelets = G(xx);

                wavelet_norm = sqrt(sum(raw_wavelets(:).^2));

                wavelets(:, k) = raw_wavelets./wavelet_norm;
            end
            
            gabor_raw = sum(repmat(AA', [size(wavelets, 1), 1]).*wavelets, 2);
            gabor_norm = sqrt(sum(gabor_raw.^2));
            gabor_coefficients = gabor_raw./gabor_norm;
        end
        
        
        
%         function [ B ] = fun_predictor_gbur( p_c, a_par, data, xx, tt )
%             B = zeros(length(xx), length(tt));
% 
%             %
%             % Convolution method
%             %
% 
%             data_box_xx = (min(xx) - (size(p_c.gabor_coeff, 1) - 1)/2):...
%                           (max(xx) + (size(p_c.gabor_coeff, 1) - 1)/2);
% 
%             data_box_xx = data_box_xx + a_par.x_max*(data_box_xx < 1);
%             data_box_xx = data_box_xx - a_par.x_max*(data_box_xx > a_par.x_max);
%             dt = ceil((size(p_c.gabor_coeff, 2) - 1)/2);
%             data_box_tt = (min(tt) - dt):(max(tt) + dt);
%             %data_box_tt = floor(data_box_tt);
% 
%             %box = abs(data(data_box_xx, data_box_tt)).^2;
%             box = data(data_box_xx, data_box_tt);
% 
%             B = conv2(box, p_c.gabor_coeff, 'valid');
%         end
    end
    
end

