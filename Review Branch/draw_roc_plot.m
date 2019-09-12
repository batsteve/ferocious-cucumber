function [outcode] = draw_roc_plot(a_par, r_par, w_bc)
%DRAW_ROC_PLOT Summary of this function goes here
%   Detailed explanation goes here

    curplot_num = a_par.next_misc_plot_num;
    if a_par.iterate_plots
        a_par.next_misc_plot_num = a_par.next_misc_plot_num + 2;
    end
    
    figure(curplot_num);
    clf;
    plot(w_bc.recall, w_bc.precision, 'LineWidth', 3);
    xlabel('recall');
    ylabel('precision');
    title('Precision-Recall ROC plot')
    set(gca, 'FontSize', 16);
    
    figure(curplot_num + 1);
    clf;
    plot(1 - w_bc.specificity, w_bc.recall, 'LineWidth', 3);
    xlabel('false positive rate');
    ylabel('true positive rate');
    title('Sensitivity-Specificity ROC plot')
    set(gca, 'FontSize', 16);
    
    drawnow();

    outcode = 1;
end

