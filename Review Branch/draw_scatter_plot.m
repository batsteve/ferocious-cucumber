function [ outcode ] = draw_scatter_plot(  a_par, r_par, h_c  )
%PLOT_SCATTER_PLOT Summary of this function goes here
%   Detailed explanation goes here

    %pre = h_c.P_hist_B;
    %ind = h_c.P_hist_A;

    
%     figure(1);
%     clf;
%     set(gcf, 'Renderer', a_par.renderer);
%     scatter(ind(:), pre(:), 1);
%     title('scatter plot of interior points')
%     xlabel('indicator')
%     ylabel('predictor')

%     h = histogram2(ind(:), pre(:), ...
%         [100, 100]);
        %[f_par.n_A, f_par.n_B]);
%     XX = h.XBinEdges(2:end);
%     YY = h.YBinEdges(2:end);
%     CC = h.Values;

    [CC, X, Y] = histcounts2(h_c.P_hist_A(:), h_c.P_hist_B(:), [48, 48]);
    [XX, YY] = meshgrid(X(2:end), Y(2:end));

%     figure(2);
%     clf;
%     set(gcf, 'Renderer', a_par.renderer);
%     title(sprintf('histogram of data points:  %s', r_par.title_tag));
%     xlabel('indicator');
%     ylabel('predictor');
    
    
    
%     figure(3);
%     clf;
%     set(gcf, 'Renderer', a_par.renderer);
%     colormap([1, 1, 1; parula(127)]);
%     pcolor(XX, YY, CC');
%     %imagesc(XX, flipud(YY), flipud(log(CC')));
%     %shading interp
%     shading flat
%     %shading faceted
%     title(sprintf('histogram:  %s', r_par.title_tag));
%     xlabel('A -- Indicator');
%     ylabel('B -- Predictor');
%     j = gca;
%     j.FontSize = 16;
%     if a_par.save_intermediate_figs
%         filename = sprintf('%sscatter_plot%s.jpg', a_par.plot_path, r_par.file_tag);
%         print('-f3', '-painters', filename,'-djpeg');
%         filename = sprintf('%sscatter_plot%s.fig', a_par.figure_path, r_par.file_tag);
%         savefig(figure(3), filename, 'compact')
%     end
    
    
    
    YYp = zeros(length(h_c.A_thresholds), length(h_c.B_thresholds));
    
    
    Ap = discretize(h_c.P_hist_A, h_c.A_thresholds);
    Bp = discretize(h_c.P_hist_B, h_c.B_thresholds);
    for k = 1:length(h_c.P_hist_N)
        YYp(Ap(k), Bp(k)) = YYp(Ap(k), Bp(k)) + h_c.P_hist_N(k);
    end
    
    curplot_num = a_par.next_scatter_plot_num;
    if a_par.iterate_plots
        a_par.next_scatter_plot_num = a_par.next_scatter_plot_num + 1;
    end
     
    if a_par.use_identical_histograms
        figure(curplot_num);
        clf;
        set(gcf, 'Renderer', a_par.renderer);
        %pcolor(h_c.A_thresholds, h_c.B_thresholds, log(YYp)');
        pcolor(h_c.A_thresholds, h_c.B_thresholds, YYp');
        shading flat
        %title(sprintf('histogram:  %s', r_par.title_tag));
        %xlabel('A -- Indicator');
        %ylabel('B -- Predictor');
        title('prediction-truth joint pdf', 'Interpreter', 'Latex')
        xlabel('Indicator ($A$)', 'Interpreter', 'Latex');
        ylabel('Predictor ($B$)', 'Interpreter', 'Latex');
        j = gca;
        j.FontSize = 14;
    else
        YYr = zeros(length(h_c.A_thresholds), length(h_c.B_thresholds));
        YYf = zeros(length(h_c.A_thresholds), 1);
        
        Ar = discretize(h_c.R_hist_A, h_c.A_thresholds);
        Br = discretize(h_c.R_hist_B, h_c.B_thresholds);
        for k = 1:length(h_c.P_hist_N)
            YYr(Ar(k), Br(k)) = YYr(Ar(k), Br(k)) + h_c.R_hist_N(k);
        end

        Af = discretize(h_c.F_hist_A, [0, h_c.A_thresholds, inf]);
        Af(Af > length(h_c.A_thresholds)) = length(h_c.A_thresholds);
        for k = 1:length(h_c.F_hist_N)
            YYf(Af(k)) = YYf(Af(k)) + h_c.F_hist_N(k);
        end
        
        figure(curplot_num);
        clf;
        set(gcf, 'Renderer', a_par.renderer);
        subplot(2, 2, 1);
        pcolor(h_c.A_thresholds, h_c.B_thresholds, log(YYp)');
        shading flat
        title(sprintf('precision hist:  %s', r_par.title_tag));
        xlabel('A -- Indicator');
        ylabel('B -- Predictor');
        j = gca;
        j.FontSize = 16;
        subplot(2, 2, 2);
        pcolor(h_c.A_thresholds, h_c.B_thresholds, log(YYr)');
        shading flat
        title(sprintf('recall hist:  %s', r_par.title_tag));
        xlabel('A -- Indicator');
        ylabel('B -- Predictor');
        j = gca;
        j.FontSize = 16;
        subplot(2, 2, 3);
        plot(h_c.A_thresholds, YYf);
        shading flat
        title(sprintf('ferocity hist:  %s', r_par.title_tag));
        xlabel('A -- Indicator');
        j = gca;
        j.FontSize = 16;
    end
    
    if a_par.save_intermediate_figs
        filename = sprintf('%spdf_hist%s.jpg', a_par.plot_path, r_par.file_tag);
        print(sprintf('-f%d', curplot_num), '-painters', filename,'-djpeg');
        filename = sprintf('%spdf_hist%s.fig', a_par.figure_path, r_par.file_tag);
        savefig(figure(curplot_num), filename, 'compact')
    end
    
%     figure(4);
%     clf;
%     set(gcf, 'Renderer', a_par.renderer);
%     subplot(2, 1, 1);
%     ha = histogram(ind(:), 100);
%     title('marginal histogram -- indicator')
%     j = gca;
%     j.FontSize = 16;
%     subplot(2, 1, 2);
%     hb = histogram(pre(:), 100);
%     title('marginal histogram -- predictor')
%     j = gca;
%     j.FontSize = 16;
    
    
    drawnow();
    
    outcode = 1;

end

