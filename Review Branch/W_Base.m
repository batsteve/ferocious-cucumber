classdef W_Base < handle
    %W_BASE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        w_stat;
        w_prc;
        w_pca;
        w_bc;
    end
    
    methods
        function [ w_base ] = W_Base()
            w_base.w_stat = 0;
            w_base.w_prc = 0;
            w_base.w_pca = 0;
            w_base.w_bc = 0;
        end
    end
    
end

