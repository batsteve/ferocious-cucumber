function [outcode] = test_plot(par)
%TEST_PLOT Summary of this function goes here
%   Detailed explanation goes here

    load(par.indicator_filename, '-mat');
    
    TT = par.t_start:par.t_step:par.t_end;
    TT = TT(2:end)';
    
    figure(11);
    clf;
    hold on
    plot(TT, e_diss_mu./max(e_diss_mu), 'Color', 'Red');
    plot(TT, abs(f_psi(:, 2, 2))./max(abs(f_psi(:, 2, 2))), 'Color', 'Blue');
    legend('loss', 'inflow');
    hold off


    outcode = 1;
end

