function w=vort(u1,u2)
global kx ky

fu1=fft2(u1);
fu2=fft2(u2);
fw=complex(0,kx).*fu2-complex(0,ky).*fu1;
w=ifft2(fw,'symmetric');