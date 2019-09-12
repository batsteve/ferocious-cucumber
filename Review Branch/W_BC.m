classdef W_BC < handle
    %W_PRC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        A_threshold;
        B_threshold; 
        
        p_true_positive;
        p_true_negative;
        p_false_positive;
        p_false_negative;
        
        r_true_positive;
        r_true_negative;
        r_false_positive;
        r_false_negative;
        
        precision;
        recall;
        specificity;
        total_error_rate;
        
        balanced_error_rate;
        F_score;
        
        best_ter;
        best_ber;
        best_f_score;
        
    end
    
    methods
        function [ w_bc ] = W_BC()
            w_bc.A_threshold = 0;
            w_bc.B_threshold = 0;

            w_bc.p_true_positive = 0;
            w_bc.p_true_negative = 0;
            w_bc.p_false_positive = 0;
            w_bc.p_false_negative = 0;
            
            w_bc.r_true_positive = 0;
            w_bc.r_true_negative = 0;
            w_bc.r_false_positive = 0;
            w_bc.r_false_negative = 0;0;

            w_bc.precision = 0;
            w_bc.recall = 0;
            w_bc.specificity = 0;
            w_bc.total_error_rate = 0;
            
            w_bc.balanced_error_rate = 0;
            w_bc.F_score = 0;
            
            w_bc.best_ter = 0;
            w_bc.best_ber = 0;
            w_bc.best_f_score = 0;
        end
    end
    
end

