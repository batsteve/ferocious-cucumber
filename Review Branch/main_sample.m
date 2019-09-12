% ~128 is baseline.  Go to 256 for crisper pictures
n_thresholds = 128;
n_interp = 129;

%calc_rule = 'from_pdf';
calc_rule = 'from_simulation';

% for 'from_simulation', determines the number of simulated samples.
% Decrease this to see the effects of sampling noise
n_sims = 1e6;

% 10 and 3 are the default (stringent!) limits
q_fac = 8;  % how far out does integration go?
t_fac = 3;  % how far out do thresholds go?

draw_plots = true;

histogram_shape = 'bimodal';
%histogram_shape = 'multi_gaussian';

% construct the model pdf

switch histogram_shape
    case 'bimodal'
        gamma = 0.2;    % weighting factor
        r = 2;        % mode separation
        
        par1_axis_label = '\gamma';
        par2_axis_label = '\rho';
        model_title = 'bimodal';
        
        pdf = @(x, y, gamma, r) 1*exp(-(x.^2 + y.^2).*r.^2) + ...
            gamma.*exp(-((x-1).^2 + (y-1).^2).*r.^2);
        
        a_min_lim_rule = @(alpha, r) -q_fac/r;
        a_max_lim_rule = @(alpha, r) 1 + q_fac/r;
        b_min_lim_rule = @(alpha, r) -q_fac/r;
        b_max_lim_rule = @(alpha, r) 1 + q_fac/r;
        
        a_thresh_rule = @(alpha, r) linspace(-t_fac/r, 1+t_fac/r, 2*n_thresholds+1);
        b_thresh_rule = @(alpha, r) linspace(-t_fac/r, 1+t_fac/r, 2*n_thresholds+1);
        
        par1 = gamma;
        par2 = r;
        
    case 'multi_gaussian'
        
        r = 0.2;        % scale factor
        theta = pi/3;   % principle axis angle
        
        par1_axis_label = '\rho';
        par2_axis_label = '\theta';
        model_title = 'multi-gaussian';
        
        s = @(r, theta) abs(r*cos(theta)) + 1./r.*sin(theta);
        t = @(r, theta) r*sin(theta) + abs(1./r.*cos(theta));
        
        pdf = @(x, y, r, theta) exp(-((x.*cos(theta) + y.*sin(theta))./r).^2 - ...
            ((x.*sin(theta) - y.*cos(theta)).*r).^2);
        
        a_min_lim_rule = @(r, theta) -q_fac*s(r, theta);
        a_max_lim_rule = @(r, theta)  q_fac*s(r, theta);
        b_min_lim_rule = @(r, theta) -q_fac*t(r, theta);
        b_max_lim_rule = @(r, theta)  q_fac*t(r, theta);
        
        a_thresh_rule = @(r, theta) linspace(-t_fac*s(r, theta),t_fac*s(r, theta), 2*n_thresholds+1);
        b_thresh_rule = @(r, theta) linspace(-t_fac*t(r, theta),t_fac*t(r, theta), 2*n_thresholds+1);
        
        par1 = r;
        par2 = theta;
        
    otherwise
        warning('Model not recognized!');
end


tic

a_min = a_min_lim_rule(par1, par2);
a_max = a_max_lim_rule(par1, par2);
b_min = b_min_lim_rule(par1, par2);
b_max = b_max_lim_rule(par1, par2);

total_weight = integral2(@(a, b) pdf(a, b, par1, par2), a_min, a_max, b_min, b_max);
pdf_norm = @(a, b, alpha, r) pdf(a, b, par1, par2)/total_weight;


%
% calculate discrete probabilities over a convenient grid
%

X = a_thresh_rule(par1, par2);
Y = b_thresh_rule(par1, par2);
[XX, YY] = meshgrid(X, Y);

A = @(a1, a2, b1, b2, alpha, r) integral2(@(x, y) pdf_norm(x, y, alpha, r), a1, a2, b1, b2);
AA = zeros(size(XX));

Xs = [a_min, a_thresh_rule(par1, par2), a_max];
Ys = [b_min, b_thresh_rule(par1, par2), b_max];
[XXs, YYs] = meshgrid(Xs, Ys);

for k = 1:(size(XXs, 1) - 1)
    for j = 1:(size(XXs, 2) - 1)
        AA(k, j) = A(Xs(k), Xs(k+1), Ys(j), Ys(j+1), par1, par2);
    end
end



switch calc_rule
    case 'from_pdf'
        % calculate analytic quantities, directly from the cumulative
        % probability functions

        BB = AA;
        
    case 'from_simulation'
        % calculate simulated quantities, from sampled cumulative
        % probability functions

        cdf = cumsum(AA(:));
        cdf = cdf./max(cdf);
        cdf = [cdf; inf];
        
        p_samples = rand(n_sims, 1);
        B_linear = histcounts(p_samples, cdf);
        
        BB = reshape(B_linear, (size(XXs, 1) - 1), (size(XXs, 2) - 1));
        BB = BB./sum(BB(:));
          
end

%
% Calculate derived quantities using cumulative probability functions
%

DD = zeros(size(XX));
EE = zeros(size(X));
FF = zeros(size(Y));

for k = 1:length(X)
    for j = 1:length(Y)
        DD(k, j) = sum(sum(BB((1+k):end, (1+j):end), 2), 1);
    end
end
for k = 1:length(X)
    EE(k) = sum(sum(BB((1+k):end, 1:end), 2), 1);
end
for j = 1:length(Y)
    FF(j) = sum(sum(BB(1:end, (1+j):end), 2), 1);
end

%
% precision, recall, and rate (abbreviated 'Fer' for 'ferocity')
%

Pre = DD./(FF);
Rec = DD./(EE');
Fer = repmat(EE', [1, size(XX, 1)]);

%
% Polish some sharp edges off (usually a result of dividing by zero
% somewhere)
%

Pre(Pre > 1) = 1;
Rec(Rec > 1) = 1;
Fer(Fer > 1) = 1;

Pre(Pre < 0) = 0;
Rec(Rec < 0) = 0;
Fer(Fer < 0) = 0;

Pre(isnan(Pre)) = 0;
Rec(isnan(Rec)) = 0;
Fer(isnan(Fer)) = 0;

%
% Interpolate onto QRS grid
%

pR = linspace(0, 1, n_interp);
pQ = linspace(0, 1, n_interp);
[pRR, pQQ] = meshgrid(pR, pQ);
II = scatteredInterpolant(Rec(:), Fer(:), Pre(:), 'natural');
pSS = II(pRR, pQQ);

pSS(isnan(pSS)) = 0;

%
% Plots!
%

if draw_plots
    
    density = pdf(XX, YY, par1, par2);
    
    figure(1);
    clf;
    set(gcf, 'Renderer', 'painters');
    mesh(XX, YY, density);
    xlabel('$A$ (indicator)', 'Interpreter', 'Latex');
    ylabel('$B$ (predictor)', 'Interpreter', 'Latex');
    zlabel('probability density', 'Interpreter', 'Latex');
    title(sprintf('pdf %s ($%s = %0.2f$, $%s= %0.1f$)', ...
        model_title, par1_axis_label, par1, par2_axis_label, par2), ...
        'Interpreter', 'Latex');
    set(gca, 'FontSize', 16);
    
    figure(2);
    clf;
    set(gcf, 'Renderer', 'painters');
    colormap([1 1 1; parula(15)]);
    pcolor(XX, YY, density);
    shading flat
    xlabel('$A$ (indicator)', 'Interpreter', 'Latex');
    ylabel('$B$ (predictor)', 'Interpreter', 'Latex');
    zlabel('probability density', 'Interpreter', 'Latex');
    title(sprintf('pdf %s ($%s = %0.2f$, $%s= %0.1f$)', ...
        model_title, par1_axis_label, par1, par2_axis_label, par2), ...
        'Interpreter', 'Latex');
    set(gca, 'FontSize', 16);
    
    
    if n_thresholds > 32
        point_size = 5;
    else
        point_size = 36;
    end
    
    figure(3);
    clf;
    set(gcf, 'Renderer', 'painters');
    scatter3(Rec(:), Fer(:), Pre(:), point_size);
    xlabel('$r$', 'Interpreter', 'Latex');
    ylabel('$q$', 'Interpreter', 'Latex');
    zlabel('$s$', 'Interpreter', 'Latex');
    xlim([0, 1]);
    ylim([0, 1]);
    zlim([0, 1]);
    view(45,30)
    title(sprintf('QRS surface: %s ($%s = %0.2f$, $%s= %0.1f$)', ...
        model_title, par1_axis_label, par1, par2_axis_label, par2), ...
        'Interpreter', 'Latex');
    set(gca, 'FontSize', 16);
    
    
    figure(4);
    clf;
    set(gcf, 'Renderer', 'painters');
    mesh(pRR, pQQ, pSS);
    xlabel('$r$', 'Interpreter', 'Latex');
    ylabel('$q$', 'Interpreter', 'Latex');
    zlabel('$s$', 'Interpreter', 'Latex');
    xlim([0, 1]);
    ylim([0, 1]);
    zlim([0, 1]);
    view(45,30)
    title(sprintf('QRS surface: %s ($%s = %0.2f$, $%s= %0.1f$)', ...
        model_title, par1_axis_label, par1, par2_axis_label, par2), ...
       'Interpreter', 'Latex');
    set(gca, 'FontSize', 16);
    
    drawnow();
end

%
% Calculate summary statistics
%

dQ = 1/(n_interp - 1);
dR = 1/(n_interp - 1);

pSS(isnan(pSS)) = 0;

center_SS = pSS(2:(end-1), 2:(end-1));
side_SS = [pSS(2:(end-1), 1); pSS(2:(end-1), end); ...
    pSS(1, 2:(end-1))';  pSS(end, 2:(end-1))'];
corner_SS = [pSS(1, 1); pSS(1, end); pSS(end, 1);  pSS(end, end)];

tot_VUS = (sum(center_SS(:)) + 1/2*sum(side_SS(:)) + 1/4*sum(corner_SS(:))) *...
    (dQ*dR);

%adjusted_auc = (1/2*pPP(1, :) + sum(pPP(2:(end-1), :), 1) + 1/2*pPP(end, :))*dF - pF;
adjusted_auc = (1/2*pSS(:, 1) + sum(pSS(:, 2:(end-1)), 2) + 1/2*pSS(:, end))*dQ - pQ';
max_adjusted_auc = max(adjusted_auc(:));

fprintf('  Took %0.2f seconds.\n', toc);

fprintf('Total volume under the surface (V):  %0.2f.\n', tot_VUS);
fprintf('Maximum adjusted area under the curve (alpha*):  %0.2f.\n', max_adjusted_auc);
