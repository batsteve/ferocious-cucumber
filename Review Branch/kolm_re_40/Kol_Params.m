classdef Kol_Params < handle
    %KOL_PARAMS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        save_vector_fields;
        data_folder;
        data_filename_base;
        indicator_filename;
        
        L1;
        L2;
        n1;
        n2;
        
        t_start;
        t_end;
        t_step;
        ode45_options;
        
        initial_conditions;
        input_file;
        
        re;
        nu;
        forcing_amplitude;
        forcing_wavenumber;
        
        plot_intermediate_vector_fields;
    end
    
    methods
        function par = Kol_Params()
            par.save_vector_fields = false;
            par.data_folder = 'data';
            
            if ~exist(par.data_folder, 'dir')
                mkdir(par.data_folder);
            end
            
            par.data_filename_base = sprintf('%s/turb_u_', par.data_folder);
            par.indicator_filename = sprintf('%s/full_data.dat', par.data_folder);
            
            par.L1=2*pi; par.n1=128;
            par.L2=2*pi; par.n2=128;
            
            par.t_start = 0;
            par.t_end = 20;
            par.t_step = 1;
            
            par.ode45_options = odeset('RelTol',1e-5,'AbsTol',1e-5);
            
            par.re = 40;
            par.nu = 1/par.re;
            
            par.forcing_amplitude = 1;
            par.forcing_wavenumber = 4;
            
            par.initial_conditions = 2;
            par.input_file = 0;
            
            par.plot_intermediate_vector_fields = false;
        end
    end
end

