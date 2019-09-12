function [u1, u2]=u0(par)
    global n1 n2 kx ky x y nu

    switch par.initial_conditions
        case 1 %read from file
            wfile=['turb_u_' sprintf('%4.4i', par.inputfile)];
            load(wfile, 'u1','u2');

        case 2 %random with gaussian decay
            k=sqrt(kx.^2+ky.^2); index=find(k==0); k(index)=1;
            k0=2;
            phi=rand(n2,n1);
            E=k.^4.*exp(-(k./k0).^2);
            E0=2e1*sqrt(2)/sum(sum(E));
            fu1=E0*(-ky./k).*sqrt(E./(2*pi*k)).*exp(complex(0,2*pi*phi));
            fu2=E0*( kx./k).*sqrt(E./(2*pi*k)).*exp(complex(0,2*pi*phi));
            fu1=fu1.*(n1*n2); fu2=fu2.*(n1*n2);
            [fu1,fu2]=dealiase(fu1,fu2);
            fu1(index)=0; fu2(index)=0;
            u1 = ifft2(fu1,'symmetric');
            u2 = ifft2(fu2,'symmetric');

        case 3 % random phase - exp decay
            k=sqrt(kx.^2+ky.^2); index=find(k==0); k(index)=1;
            k0=1;
            phi=rand(n2,n1);
            E=exp(-k./k0);
            E0=-sqrt(2);
            fu1=E0*(-ky./k).*sqrt(E./(2*pi*k)).*exp(complex(0,2*pi*phi));
            fu2=E0*( kx./k).*sqrt(E./(2*pi*k)).*exp(complex(0,2*pi*phi));

            % set I=I_0
            I0=.1;
            fu1(5,1)=-complex(0,I0);
            fu1(125,1) = conj(fu1(5,1));

            fu1=fu1*(n1*n2); fu2=fu2*(n1*n2);
            [fu1,fu2]=dealiase(fu1,fu2);
            fu1(index)=0; fu2(index)=0;
            % Fourier to Phys.
            u1 = ifft2(fu1,'symmetric');
            u2 = ifft2(fu2,'symmetric');

        case 4 % exponential decay - uniformly random phase - gaussian distributed amp
            k=sqrt(kx.^2+ky.^2); index=find(k==0); k(index)=1;
            k0=2;
            rng('shuffle');
            phi=rand(n2,n1);
            Phi=normrnd(0*k,exp(-k./k0));
            E0=sqrt(2);
            fu1=E0*(-ky./k).*sqrt(Phi.^2./(2*pi*k)).*exp(complex(0,2*pi*phi));
            fu2=E0*( kx./k).*sqrt(Phi.^2./(2*pi*k)).*exp(complex(0,2*pi*phi));

            fu1(index)=0; fu2(index)=0;

            % set I=I_0
            I0=.1;
            fu1(5,1)=-complex(0,I0);
            fu1(125,1) = conj(fu1(5,1));

            fu1=fu1*(n1*n2); fu2=fu2*(n1*n2);
            [fu1,fu2]=dealiase(fu1,fu2);

            % Fourier to Phys.
            u1 = ifft2(fu1,'symmetric');
            u2 = ifft2(fu2,'symmetric');

        case 5
            u1 = (1/(16*nu))*sin(4*y);
            u2 = 0*u1;


        otherwise 
            warning('%s not recognized!\n', par.initial_conditions);
    end

    wfile=[par.data_filename_base, sprintf('%4.4i',0)];
    if ~exist(wfile, 'file') && par.save_vector_fields
        save(wfile, 'u1','u2');
    end
end



% function [u1,u2]=gvortex(xx,yy,x0,y0,delta,gamma0)
% r2=(xx-x0).^2+(yy-y0).^2;
% index=find(r2==0);
% u1 = gamma0*(1-exp(-r2/delta^2))./(2*pi*r2).*(y0-yy);
% u1(index)=gamma0/(2*pi);
% u2 = gamma0*(1-exp(-r2/delta^2))./(2*pi*r2).*(xx-x0);
% u2(index)=gamma0/(2*pi);
% end
