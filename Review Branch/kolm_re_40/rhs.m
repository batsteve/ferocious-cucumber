function [rhs1, rhs2]=rhs(u1,u2, par)
global kx ky nu

fu1=fft2(u1); fu2=fft2(u2);

% fu1x = complex(0,kx).*fu1; u1x = ifft2(fu1x,'symmetric');
% fu1y = complex(0,ky).*fu1; u1y = ifft2(fu1y,'symmetric');
% fu2x = complex(0,kx).*fu2; u2x = ifft2(fu2x,'symmetric');
% fu2y = complex(0,ky).*fu2; u2y = ifft2(fu2y,'symmetric');
% rhs1 = -(u1.*u1x+u2.*u1y); rhs2 = -(u1.*u2x+u2.*u2y);

ksq = kx.^2+ky.^2;

w=vort(u1,u2);
rhs1 =  u2.*w; rhs2 = -u1.*w;

frhs1 = fft2(rhs1); frhs2=fft2(rhs2);

% add force + drag
[frc1, frc2]=force(par);
frhs1=frhs1+frc1;
frhs2=frhs2+frc2;

% dissipation
frhs1 = frhs1- nu*ksq.*fu1; 
frhs2 = frhs2- nu*ksq.*fu2;

% dealiase + project to div free fields
[frhs1, frhs2] = dealiase(frhs1,frhs2);
[frhs1, frhs2] = divfree(frhs1,frhs2);

rhs1 = ifft2(frhs1,'symmetric');
rhs2 = ifft2(frhs2,'symmetric');