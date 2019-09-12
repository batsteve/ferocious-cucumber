classdef Run_Parameters < handle
    %RUN_PARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        analysis_state;
        run_number;
        
        file_tag;
        title_tag;
        
        % kolm paramaters
        knx;
        kny;
        kx;
        ky;
        
        % mmt parameters
        gabor_L_real;
        gabor_L_steps;
        kernel;
        
        % arbitrary optimization parameters
        coeffs;
        pred_fun;
        
        % test parameters
        alpha;
        delta_r;
        
        % danger flags
        dangerous_prc;
    end
    
    methods
        function [ r_par ] = Run_Parameters()
            r_par.analysis_state = -1;
            r_par.run_number = 0;
            
            r_par.file_tag = '';
            r_par.title_tag = '';
            
            r_par.knx = 0;
            r_par.kny = 0;
            r_par.kx = 0;
            r_par.ky = 0;
            
            r_par.gabor_L_real = 0;
            r_par.gabor_L_steps = 0;
            r_par.kernel = '';
            
            r_par.coeffs = 0;
            r_par.pred_fun = [];
            
            r_par.alpha = 0;
            r_par.delta_r = 4;
            
            r_par.dangerous_prc = false;
        end
    end
    
end

