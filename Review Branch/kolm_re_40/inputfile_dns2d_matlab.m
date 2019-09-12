clc;
clear;
global n1 n2 L1 L2 kx ky x y nu
L1=2*pi; n1=128;
L2=2*pi; n2=128;

[x,y,kx,ky]=gvars(n1,n2,L1,L2);

tinit =0;
tend  =20;
dwrite=1;
ic=2;
Re=40; nu=1/Re;

[u1, u2]=u0(ic,0);
[u1,u2] = dns2d(u1,u2,tinit,tend,dwrite);

plot_turb(tinit,tend,dwrite)