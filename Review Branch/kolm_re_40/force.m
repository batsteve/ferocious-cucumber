function [frc1,frc2]=force(par)
global y

a1 = par.forcing_amplitude*sin(par.forcing_wavenumber*y);
a2 = 0*a1;
frc1 = fft2(a1);
frc2 = fft2(a2);