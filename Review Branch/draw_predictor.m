set(groot, 'DefaultFigureRenderer', 'painters');
set(0, 'DefaultFigureRenderer', 'painters');

%
%
%
coeffs = mean(squeeze(best_x(4, 2, :, :)), 1);
%
%
%

a_par = Analysis_Parameters();
addpath('../../MMT_NL/main_version');
a_par.datapath = '../../../Data/MMT/July';
a_par.gap_T_real = 0.015;
a_par.analysis_name = sprintf('_t_%0.2f', a_par.gap_T_real);
a_par.verbosity = 2;
filename = sprintf('%s/sim_data/sim_%d.dat', a_par.datapath, 1);
fprintf('Loading file %s in order to read data parameters.\n', filename);

load(filename, 'm_par', '-mat');

a_par.drop_T_steps = 1000;
a_par.gap_T_steps = ceil(a_par.gap_T_real./m_par.dt);
a_par.gap_T_real = a_par.gap_T_steps*m_par.dt;
a_par.analysis_name = sprintf('_t_%0.4f', a_par.gap_T_real);

a_par.evaluation_mode = 'test';
a_par.nonlinear_constraint_type = 'none';
%a_par.testing_file_set = 21:28;
a_par.testing_file_set = 23:30;
a_par.hypothesis_class_size = 2;

r_par = Run_Parameters();
a_par.cur_run = a_par.cur_run + 1;
r_par.run_number = a_par.cur_run;
r_par.coeffs = coeffs;
r_par.file_tag = sprintf('_n=%d', r_par.run_number);
r_par.title_tag = sprintf('(n=%d)', r_par.run_number);
r_par.pred_fun = Predictor_Container(a_par, coeffs);

[ AA, LL ] = parametersToWavelets( a_par, coeffs );
h_c = Histogram_Container(a_par, r_par);

[CC, X, Y] = histcounts2(h_c.P_hist_A(:), h_c.P_hist_B(:), [48, 48]);
[XX, YY] = meshgrid(X(2:end), Y(2:end));

YYp = zeros(length(h_c.A_thresholds), length(h_c.B_thresholds));


Ap = discretize(h_c.P_hist_A, h_c.A_thresholds);
Bp = discretize(h_c.P_hist_B, h_c.B_thresholds);
for k = 1:length(h_c.P_hist_N)
    YYp(Ap(k), Bp(k)) = YYp(Ap(k), Bp(k)) + h_c.P_hist_N(k);
end

[x, y] = ndgrid(linspace(-3, 3, 9), linspace(-3, 3, 9));
sm = exp(-x.^2 + -y.^2);
sm = sm./(sum(sm(:)));
YYq = conv2(YYp, sm, 'same');

figure(11);
clf;
set(gcf, 'Renderer', a_par.renderer);
imagesc(h_c.A_thresholds, h_c.B_thresholds, YYp');
%shading flat
title('prediction-truth joint pdf');
xlabel('Truth');
ylabel('Prediction');
xlim([0.5 3.2])
j = gca;
j.FontSize = 16;
j.YDir = 'normal';
j.CLim = [0, 7];


a_par.save_intermediate_figs = 0;
w_base = W_Base();
compute_prc_quantities( h_c, a_par, w_base );
draw_prc_plot(a_par, r_par, w_base.w_prc);