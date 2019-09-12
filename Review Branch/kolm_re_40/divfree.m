function [f1 f2]=divfree(f1,f2);
global kx ky

tmp1=f1; tmp2=f2;

ksq = kx.^2+ky.^2;
index = find(ksq==0);
ksq(index)=1;

f1 = (ky.^2.*tmp1 - kx.*ky.*tmp2)./ksq;
f2 = (-kx.*ky.*tmp1 + kx.^2.*tmp2)./ksq;