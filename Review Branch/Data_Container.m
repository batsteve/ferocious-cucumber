classdef Data_Container < handle
    %DATA_CONTAINER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % mmt
        raw_sim;
        
        % kolm
        
        raw_ind;
        raw_pre;
    end
    
    methods
        function [ data_c ] = Data_Container()
            data_c.raw_sim = 0;
            
            data_c.raw_ind = 0;
            data_c.raw_pre = 0;
        end
    end
    
end

