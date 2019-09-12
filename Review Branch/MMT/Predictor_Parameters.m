classdef Predictor_Parameters < handle
    %MMT_PARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        n_x;
        
        t_start;
        t_end;
        
        gabor_L;
        
        predictor_kernel;
        gabor_coefficients;
    end
    
    methods
        function [ ab_par ] = Predictor_Parameters( m_par )
            ab_par.n_x = m_par.n_x;
                        
            ab_par.t_start = 1000;
            ab_par.t_end = m_par.n_saved - 100;
            
            ab_par.gabor_L = 1;
            
            ab_par.predictor_kernel = 'square_exp';
            %ab_par.predictor_kernel = 'sinc';
            ab_par.gabor_coefficients = 0;
        end
    end
    
end

