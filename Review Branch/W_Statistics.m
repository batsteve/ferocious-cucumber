classdef W_Statistics < handle
    %W_STATISTICS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        full_counts;
        full_A_mean;
        full_B_mean;
        full_A_var;
        full_B_var;
    end
    
    methods
        function [ w_stat ] = W_Statistics()
            w_stat.full_counts = 0;
            w_stat.full_A_mean = 0;
            w_stat.full_B_mean = 0;
            w_stat.full_A_var = 0;
            w_stat.full_B_var = 0;
        end
    end
    
end

