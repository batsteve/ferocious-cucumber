classdef Hist_Parameters
    %MMT_PARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        n_x;
        n_bins;
        
        t_start;
        t_end;
        %x_decimate;
        %t_decimate;
    end
    
    methods
        function [ h_par ] = Hist_Parameters( ab_par )
            %m_par.n_sims = 1;
            
            % It uses Cox and Matthews's ETDRK4 exponential integrator with a fixed
            % time step.
            h_par.n_x = ab_par.n_x;
            h_par.n_bins = 4000;
            
            h_par.t_start = ab_par.t_start;
            h_par.t_end = ab_par.t_end;
            
            %m_par.x_decimate = f_par.x_decimate;
            %m_par.t_decimate = f_par.t_decimate;
        end
    end
    
end

