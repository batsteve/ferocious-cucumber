function [u1, u2]=dns2d(u1,u2,par)
global L1 L2 n1 n2 x y

%
% build data structures
%

u1_vec = reshape(u1,n1*n2,1);
u2_vec = reshape(u2,n1*n2,1);

kk_full = -7:7;
kk_half = 0:7;
[xx, yy] = ndgrid(2*pi*(1:n1)./n1, 2*pi*(1:n2)/n2);

TT = par.t_start:par.t_step:(par.t_end-par.t_step);
e_diss_mu = zeros(length(TT), 1);
e_kinetic_mu = zeros(length(TT), 1);
e_inflow_mu = zeros(length(TT), 1);
vort_mu = zeros(length(TT), 1);
q_mu = zeros(length(TT), 1);

f_u1 = zeros([length(TT), length(kk_full), length(kk_half)]);
f_u2 = zeros([length(TT), length(kk_full), length(kk_half)]);
f_u12 = zeros([length(TT), length(kk_full), length(kk_half)]);
f_psi = zeros([length(TT), length(kk_full), length(kk_half)]);

% strain tensor eigenvalues?
% McIlhaney beta? (or alpha, or gamma)
% local hessian quantities?
% FTLE?

% time derivative
% moving average



%
% set initial conditions
%

%A = L1*L2/(n1*n2);
[kF1, kF2] = force(par);
F1 = ifft2(kF1);
F2 = ifft2(kF2);

tic
fprintf('Beginning integration:  initial transients end at t=%0.2f.\n', ....
         par.t_start);
 
%
% integrate out the transients
%

integrand = @(t, y) oderhs(t, y, par);
[~, F]=ode45(integrand,...
                [0, par.t_start],...
                [u1_vec;u2_vec], par.ode45_options );
u1 = reshape(F(end,1:n1*n2),         n2, n1);
u2 = reshape(F(end,n1*n2+1:2*n1*n2), n2, n1);

t = 1;
calc_derived_quantities();

u1_vec = F(end,1:n1*n2)';
u2_vec = F(end,n1*n2+1:2*n1*n2)';


fprintf('Integration reached the end of the initial transient zone.\n');

for t = 2:length(TT)
    fprintf('Reached step %d (time %0.2f) after %0.2f seconds.\n', ...
        t, TT(t), toc);
    
    [~, F]=ode45(integrand,...
                [TT(t-1) TT(t-1)+par.t_step/2 TT(t)],...
                [u1_vec;u2_vec], par.ode45_options );
    u1 = reshape(F(3,1:n1*n2),         n2, n1);
    u2 = reshape(F(3,n1*n2+1:2*n1*n2), n2, n1);

    % save the vel field
    if par.save_vector_fields
        j=round(TT(t)/par.t_step)+1;
        wfile=[par.data_filename_base, sprintf('%4.4i',j)];
        save(wfile, 'u1','u2');
    end
    
    calc_derived_quantities();
    if par.plot_intermediate_vector_fields
        plot_derived_quantities
    end
    
    u1_vec = F(3,1:n1*n2)';
    u2_vec = F(3,n1*n2+1:2*n1*n2)';
end

fprintf('Integration over after %0.2f seconds!\n', toc);

u1 = reshape(F(3,1:n1*n2),         n2, n1);
u2 = reshape(F(3,n1*n2+1:2*n1*n2), n2, n1);

save(par.indicator_filename, 'par', 'TT', 'e_diss_mu', 'f_u1',...
     'f_u2', 'f_u12', 'f_psi', 'e_inflow_mu', ...
     'e_kinetic_mu', 'vort_mu', 'q_mu', 'u1', 'u2');

 
 





function [] = calc_derived_quantities()
    local_e_kinetic = (u1.^2 + u2.^2);
    e_kinetic_mu(t) = mean(local_e_kinetic(:));
%     e_kinetic_var(t) = var(local_e_kinetic(:));
%     
    dudx = -1/(2*L1/n1)*([u1(end, :); u1(1:end-1,:)] - [u1(2:end, :); u1(1, :)]);
    dudy = -1/(2*L2/n2)*([u1(:, end), u1(:, 1:end-1)] -[u1(:, 2:end), u1(:, 1)]);
    dvdx = -1/(2*L1/n1)*([u2(end, :); u2(1:end-1,:)] - [u2(2:end, :); u2(1, :)]);
    dvdy = -1/(2*L2/n2)*([u2(:, end), u2(:, 1:end-1)] -[u2(:, 2:end), u2(:, 1)]);
   
    local_e_diss = par.nu*(dudx.^2 + dudy.^2 + dvdx.^2 + dvdy.^2);
    e_diss_mu(t) = mean(local_e_diss(:));
%     e_diss_var(t) = var(local_e_diss(:));
%     
    local_e_influx = (u1.*F1 + u2.*F2);
    e_inflow_mu(t) = mean(local_e_influx(:));
%     e_inflow_var(t) = var(local_e_influx(:));
%     
    local_vort = vort(u1,u2);
    vort_mu(t) = mean(local_vort(:));
%     vort_var(t) = var(local_vort(:));
%     
%     % partially unrolled tensor
%     % Frobenius norm
    local_S = 1/2*([dudx(:), dudy(:), dvdx(:), dvdy(:)] + ...
                   [dudx(:), dvdx(:), dudy(:), dvdy(:)]);
    local_S_norm = sum(abs(local_S).^2, 2);
    local_S_norm = reshape(local_S_norm, n1, n2);
    local_W = 1/2*([dudx(:), dudy(:), dvdx(:), dvdy(:)] - ...
                   [dudx(:), dvdx(:), dudy(:), dvdy(:)]);
    local_W_norm = sum(abs(local_W).^2, 2);           
    local_W_norm = reshape(local_W_norm, n1, n2);
    
    local_q = 1/2*(local_S_norm - local_W_norm);
    q_mu(t) = mean(local_q(:));
%     q_var(t) = var(local_q(:));

    psi = velocity2stream(u1, u2);

    for kn1 = 1:length(kk_full)
        for kn2 = 1:length(kk_half)
            k1 = kk_full(kn1);
            k2 = kk_half(kn2);
            
            e_fac = exp(-1i*(k1.*xx(:) + k2*yy(:)));
            f_u1(t, kn1, kn2) = sum(u1(:).*e_fac(:));
            f_u2(t, kn1, kn2) = sum(u2(:).*e_fac(:));
            f_u12(t, kn1, kn2) = sum((u1(:) + 1i*u2(:)).* e_fac(:));
            f_psi(t, kn1, kn2) = sum(psi(:).*e_fac(:));
        end
    end
end

function [] = plot_derived_quantities()
    figure(11);
    clf;
    contourf(x,y,local_e_kinetic);
    axis equal tight; colorbar
    title(sprintf('local kinetic energy at t=%0.2f', TT(t)));
    
    figure(12);
    clf;
    contourf(x,y,local_e_diss);
    axis equal tight; colorbar
    title(sprintf('local energy dissipation at t=%0.2f', TT(t)));
    
    figure(13);
    clf;
    contourf(x,y,local_e_influx);
    axis equal tight; colorbar
    title(sprintf('local energy inflow at t=%0.2f', TT(t)));
    
    figure(14);
    clf;
    contourf(x,y,local_vort);
    axis equal tight; colorbar
    title(sprintf('local vorticity at t=%0.2f', TT(t)));
    
    figure(15);
    clf;
    contourf(x,y,local_q);
    axis equal tight; colorbar
    title(sprintf('local Okubu-Weiss at t=%0.2f', TT(t)));
    
    drawnow();
end

end

