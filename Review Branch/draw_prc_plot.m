function [ outcode ] = draw_prc_plot( a_par, r_par, w_prc )
%DRAW_PRC_PLOT Summary of this function goes here
%   Detailed explanation goes here

%     figure(11)
%     clf;
%     set(gcf, 'Renderer', a_par.renderer);
%     scatter3(w_prc.recall(:), w_prc.ferocity(:), w_prc.precision(:));
%     xlabel('recall');
%     zlabel('precision');
%     ylabel('ferocity');
%     title(sprintf('3D PRF scatter plot:  %s', r_par.title_tag))
%     xlim([0, 1]);
%     ylim([0, 1]);
%     zlim([0, 1]);
%     view(45,30)
%     j = gca;
%     j.FontSize = 16;
    
    
    R = linspace(min(w_prc.recall(:)), max(w_prc.recall(:)), a_par.auc_plot_points);
    F = linspace(min(w_prc.ferocity(:)), max(w_prc.ferocity(:)), a_par.auc_plot_points);
    [RR, FF] = meshgrid(R, F);
    data_R = w_prc.recall(:);
    data_F = w_prc.ferocity(:);
    data_P = w_prc.precision(:);
    II = scatteredInterpolant(data_R, data_F, data_P);
    PP_fero = II(RR, FF);
    
    
    
    curplot_num = a_par.next_prc_plot_num;
    if a_par.iterate_plots
        a_par.next_prc_plot_num = a_par.next_prc_plot_num + 1;
    end
    
    figure(curplot_num);
    clf;
    set(gcf, 'Renderer', a_par.renderer);
    mesh(RR, FF, PP_fero);
    xlabel('recall');
    zlabel('precision');
    ylabel('ferocity');
    xlim([0, 1]);
    ylim([0, 1]);
    zlim([0, 1]);
    view(45,30)
    title(sprintf('3D PRF mesh plot:  %s', r_par.title_tag))
    j = gca;
    j.FontSize = 16;
    if a_par.save_intermediate_figs
        filename = sprintf('%sPRF_mesh_plot%s.jpg', a_par.plot_path, r_par.file_tag);
        print(sprintf('-f%d', curplot_num), '-painters', filename,'-djpeg');
        filename = sprintf('%sPRF_mesh_plot%s.fig', a_par.figure_path, r_par.file_tag);
        savefig(figure(curplot_num), filename, 'compact')
    end
    
%     figure(13)
%     clf;
%     set(gcf, 'Renderer', a_par.renderer);
%     pcolor(RR, FF, PP_fero)
%     xlabel('recall');
%     ylabel('ferocity');
%     xlim([0, 1]);
%     ylim([0, 1]);
%     caxis([0, 1]);
%     colorbar();
%     title(sprintf('PRF plot:  %s', r_par.title_tag));
%     shading flat
%     j = gca;
%     j.FontSize = 16;
%     if a_par.save_intermediate_figs
%         filename = sprintf('%sPRF_pcolor_plot%s.jpg', a_par.plot_path, r_par.file_tag);
%         print('-f13', '-painters', filename,'-djpeg');
%         filename = sprintf('%sPRF_pcolor_plot%s.fig', a_par.figure_path, r_par.file_tag);
%         savefig(figure(13), filename, 'compact')
%     end
% 
%     figure(14);
%     clf;
%     set(gcf, 'Renderer', a_par.renderer);
%     contour(RR, FF, PP_fero, 20, 'LineWidth', 3);
%     xlabel('recall');
%     ylabel('ferocity');
%     xlim([0, 1]);
%     ylim([0, 1]);
%     caxis([0, 1]);
%     colorbar();
%     title(sprintf('PRF--Precision Contours:  %s', r_par.title_tag));
%     shading flat
%     j = gca;
%     j.FontSize = 16;
%     if a_par.save_intermediate_figs
%         filename = sprintf('%sRF_contour_plot%s.jpg', a_par.plot_path, r_par.file_tag);
%         print('-f14', '-painters', filename,'-djpeg');
%         filename = sprintf('%sRF_contour_plot%s.fig', a_par.figure_path, r_par.file_tag);
%         savefig(figure(14), filename, 'compact')
%     end
%     
%     
%     %CC = [0, 0.05 0.1, 0.15 0.2, 0.25, .3, 0.35 0.5, 0.9];
%     CC = [0, 0.01, 0.05, 0.1, 0.15, 0.25, 0.5, 0.9];
%     figure(15);
%     clf;
%     set(gcf, 'Renderer', a_par.renderer);
%     contour(RR, PP_fero, FF, CC, 'LineWidth', 3);
%     xlabel('recall');
%     ylabel('precision');
%     xlim([0, 1]);
%     ylim([0, 1]);
%     caxis([0, 1]);
%     colorbar();
%     title(sprintf('PRF--Ferocity Contours:  %s', r_par.title_tag));
%     shading flat
%     j = gca;
%     j.FontSize = 16;
%     if a_par.save_intermediate_figs
%         filename = sprintf('%sPR_contour_plot%s.jpg', a_par.plot_path, r_par.file_tag);
%         print('-f15', '-painters', filename,'-djpeg');
%         filename = sprintf('%sPR_contour_plot%s.fig', a_par.figure_path, r_par.file_tag);
%         savefig(figure(15), filename, 'compact')
%     end
    
    
%     R = linspace(min(w_prc.recall(:)), max(w_prc.recall(:)), ...
%         a_par.auc_plot_points);
%     %logD = linspace(min(log(w_prc.dangerosity(:))), max(log(w_prc.dangerosity(:))), ...
%     %    a_par.auc_integration_points);
%     logD = linspace(-2, 0, a_par.auc_plot_points);
%     
%     [RR, logDD] = meshgrid(R, logD);
%     data_R = w_prc.recall(:);
%     data_logD = log(w_prc.dangerosity(:));
%     data_P = w_prc.precision(:);
%     II = scatteredInterpolant(data_R, data_logD, data_P);
%     PP_danger = II(RR, logDD);
%     
%     
%     
%     figure(16)
%     clf;
%     set(gcf, 'Renderer', a_par.renderer);
%     scatter3(w_prc.recall(:), w_prc.dangerosity(:), w_prc.precision(:));
%     xlabel('recall');
%     zlabel('precision');
%     ylabel('danger');
%     title(sprintf('3D PRD scatter plot:  %s', r_par.title_tag))
%     %xlim([0, 1]);
%     %ylim([0, 1]);
%     %zlim([0, 1]);
%     view(45,30)
%     j = gca;
%     j.FontSize = 16;
%     
%     
%     
%     figure(17);
%     clf;
%     set(gcf, 'Renderer', a_par.renderer);
%     mesh(RR, exp(logDD), PP_danger);
%     xlabel('recall');
%     zlabel('precision');
%     ylabel('danger');
%     %xlim([0, 1]);
%     %ylim([1e-2, 1]);
%     %zlim([0, 1]);
%     view(45,30)
%     title(sprintf('3D PRD mesh plot:  %s', r_par.title_tag))
%     j = gca;
%     j.FontSize = 16;
%     j.YScale = 'log';
    
    
    drawnow();

    outcode = 1;
    
end

