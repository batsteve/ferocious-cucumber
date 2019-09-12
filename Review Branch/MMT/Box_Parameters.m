classdef Box_Parameters < handle
    %MMT_PARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        n_x;
        
        t_start;
        t_end;
        
        big_delta_x;
        big_delta_t_min;
        big_delta_t_max;
        
        box_L_frac;
    end
    
    methods
        function [ b_par ] = Box_Parameters( ab_par, dt )            
            b_par.n_x = ab_par.n_x;
            
            b_par.t_start = ab_par.t_start;
            b_par.t_end = ab_par.t_end;
            
            b_par.box_L_frac = 10;
            
            b_par.big_delta_x = ceil(ab_par.gabor_L/b_par.box_L_frac);
            b_par.big_delta_t_min = abs(dt);
            b_par.big_delta_t_max = abs(dt);
        end
    end
    
end

