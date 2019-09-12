function [ YY ] = simulate_MMT_NL_data( m_par )
%SIMULATE_MMT_NL_DATA Summary of this function goes here
%   Detailed explanation goes here

    
    x = linspace(0,1, m_par.n_x+1)';
    x = x(1:m_par.n_x);
    k1=1:1:15;

    k = ([0:m_par.n_x/2 -(m_par.n_x/2-1):-1]')*(2*pi/m_par.LX);
    F = 0*0.0164*sin(4*pi*x);
    ind = find(abs(k)>500*2*pi/m_par.LX);
    D = 0*k;
    % D([2 n]) = -1; % lowest mode damping
    D(ind) = -1*(abs(k(ind))-500*2*pi/m_par.LX).^2; % Selective Laplacian high mode damping.
    clear ind

    % L = -1i*(abs(k)).^2; standard NLS
    L = -1i*sqrt(abs(k))+D;

    options = struct('lambda',m_par.lambda,'F',F,'deterministic',1==1);

    % Compute the phi functions.
    fprintf('Computing phi functions.\n')
    tic;
    beta = m_par.dt*L;
    phi00 = exp(beta);
    [phi01, phi02, phi03] = phipade(beta,3);
    phi01 =spdiags(phi01);phi02=spdiags(phi02);phi03=spdiags(phi03);
    phi10 = exp(beta/2);
    phi11 = phipade(beta/2,1);
    phi11=spdiags(phi11);
    phis = [phi00 phi01 phi02 phi03 phi10 phi11];
    clear phi0* phi1*
    fprintf('Finished phi functions after %0.2f seconds.\n', toc)
    
    fprintf('Starting simulation.\n');
    
    % Build some random initial conditions
    
    y0=0*x;
    
    for j=1:length(k1)
        y0 = y0+...
            exp(2*k1(j)*pi*1i*x+1i*rand*2*pi)+...
            exp(-2*k1(j)*pi*1i*x+1i*rand*2*pi);
    end
    y0=(y0/15+.1);
    
    %1/sqrt(20)*cos((x-0.5)*LX).*(sqrt(2)-cos((x-0.5)*LX)).^-1;
    %exp(1i*2*pi*x).*sech(100*(x-.5)); % initial condition
    
    % remove highest mode to symmetrize
    y0 = ifft(([ones(1,m_par.n_x/2) 0 ones(1,m_par.n_x/2-1)]').*fft(y0));
    
    % main loop for ETD4RK
    fprintf('Beginning integration \n')
    y=fft(y0);
    % Y = zeros(1000,n);
    tic;
    in=0; % ugh
    
    YY = zeros(m_par.n_saved, length(y));
    %yyf = zeros(n_saved, length(y));
    TT = zeros(m_par.n_saved, 1);
    
    % Simulate!
    
    for ii=1:2*m_par.n_t_max
        N0 = m_par.dt*MMT_NL(y,options);
        N1 = m_par.dt*MMT_NL(phis(:,5).*y+phis(:,6).*N0/2,options);
        N2 = m_par.dt*MMT_NL(phis(:,5).*y+phis(:,6).*N1/2,options);
        N3 = m_par.dt*MMT_NL(phis(:,1).*y+phis(:,6).*(phis(:,5)-1).*N0/2+phis(:,6).*N2,options);
        y = phis(:,1).*y+(phis*[0;1;-3;4;0;0]).*N0+(phis*[0;0;2;-4;0;0]).*N1+...
            (phis*[0;0;2;-4;0;0]).*N2+(phis*[0;0;-1;4;0;0]).*N3;
        
        if (ii>1) && (mod(ii,m_par.save_rate)==1) % only save every [save_rate] frame
            in=in+1;
            if (mod(in,1000)==0)
                fprintf('Saving step %d out of %d (%0.2f seconds).\n', in, m_par.n_saved, toc)
            end
            YY(in,:)=(ifft(y.'));
            %yyf(in,:)=0.5*(y(:)+fliplr(y(:)));
            TT(in)=ii*m_par.dt;
        end
    end
    
    YY = YY';
    fprintf('Integration finished after %02.f seconds.\n', toc)
        
end

