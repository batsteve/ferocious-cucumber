classdef MMT_Parameters < handle
    %MMT_PARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %n_sims;
        n_x;
        LX;
        Tf;
        dt;
        n_t_max;
        lambda;
        save_rate;
        n_saved;
        
        t_start;
        t_end;
        %x_decimate;
        %t_decimate;
    end
    
    methods
        function [ m_par ] = MMT_Parameters( )
            %m_par.n_sims = 1;
            
            % It uses Cox and Matthews's ETDRK4 exponential integrator with a fixed
            % time step.
            m_par.n_x = 8192;
            m_par.LX = 2*pi;
            m_par.Tf = 150;
            m_par.dt = 0.001; % dt=0.00005 for lambda=400; dt=0.001 for lambda=13;
            m_par.n_t_max = floor(m_par.Tf/m_par.dt);
            m_par.lambda = -4.0; % 1 => defocusing, -1 => focusing
            
            m_par.save_rate = 50;
            m_par.n_saved = floor(2*m_par.n_t_max / m_par.save_rate);
            
            m_par.t_start = 1000;
            m_par.t_end = m_par.n_saved - 100;
            %m_par.x_decimate = f_par.x_decimate;
            %m_par.t_decimate = f_par.t_decimate;
        end
    end
    
end

