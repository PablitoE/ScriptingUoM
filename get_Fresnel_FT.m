function FT_Fresnel = get_Fresnel_FT(dx,dy,Nc,Nr,wavel,z)
k = 2*pi/wavel;

x = - Nc/2:Nc/2-1;  % xs and ys
y = - Nr/2:Nr/2-1;
[xx,yy] = meshgrid(x*dx,y*dy);
clear('x', 'y');

% Not Fourier Transformed
FT_Fresnel = exp(1j*k*z)/(1j*wavel*z) .* exp(1j*pi/wavel/z * (xx^2+yy^2));

FT_Fresnel = fftshift(fft2(FT_Fresnel));