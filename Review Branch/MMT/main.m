n_sims = 30;
DT = [15, 30];
kernels = {'sinc', 'exp_squared'};


for k = 1:n_sims
    filename = sprintf('sim_%d.dat', k);
    YY = simulate_MMT_NL_data(m_par);
    save(filename, 'm_par', 'YY');
end
