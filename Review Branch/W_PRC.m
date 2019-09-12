classdef W_PRC < handle
    %W_PRC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        A_threshold;
        B_threshold;
        
        predictor_polarity;
        
        p_true_positive;
        p_true_negative;
        p_false_positive;
        p_false_negative;
        
        r_true_positive;
        r_true_negative;
        r_false_positive;
        r_false_negative;
        
        f_positive;
        f_negative;
        
        precision;
        recall;
        ferocity;
        specificity;
        pre_ferocity;
        dangerosity;
        
        vus;
        alpha_star;
    end
    
    methods
        function [ w_prc ] = W_PRC()
            w_prc.A_threshold = 0;
            w_prc.B_threshold = 0;
            
            w_prc.predictor_polarity = 'positive';

            w_prc.p_true_positive = 0;
            w_prc.p_true_negative = 0;
            w_prc.p_false_positive = 0;
            w_prc.p_false_negative = 0;
            
            w_prc.r_true_positive = 0;
            w_prc.r_true_negative = 0;
            w_prc.r_false_positive = 0;
            w_prc.r_false_negative = 0;
            
            w_prc.f_positive = 0;
            w_prc.f_negative = 0;

            w_prc.precision = 0;
            w_prc.recall = 0;
            w_prc.ferocity = 0;
            w_prc.specificity = 0;
            w_prc.pre_ferocity = 0;
            w_prc.dangerosity = 0;
            
            %w_prc.auc_ferocity = 0;
            %w_prc.auc_dangerosity = 0;
            
            w_prc.vus = 0;
            w_prc.alpha_star = 0;
        end
    end
    
end

