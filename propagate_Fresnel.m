function [fieldOut,xx2yy2] = propagate_Fresnel(dx,dy,wavel,z,FT_Holo,padding,xx2yy2,verbose)
if ~exist('verbose','var')
   verbose = true;
end
k = 2*pi/wavel;
[Nr, Nc] = size(FT_Holo);
dx = single(dx); dy = single(dy);
Nr = single(Nr); Nc = single(Nc);
if verbose, fprintf('Prop Fresnel: '), end
if isempty(xx2yy2)
   x = - Nc/2:Nc/2-1;  % xs and ys
   y = - Nr/2:Nr/2-1;
   [xx,yy] = meshgrid(x*dx,y*dy);
   clear('x', 'y');
   xx2yy2 = fftshift(1j*pi/wavel*(xx.^2+yy.^2));
   xx2yy2 = single(xx2yy2);
   clear('xx', 'yy');
end
if verbose, fprintf('exp term ready... '), end

% Not Fourier Transformed
FT_Fresnel = exp(1j*k*z)/(1j*wavel*z) * exp(xx2yy2 / z);
if verbose, fprintf('kernel ready... '), end

FT_Fresnel = fft2(FT_Fresnel);
FT_Holo = ifftshift(FT_Holo) .* FT_Fresnel;
clear('FT_Fresnel')
if verbose, fprintf('mult ready... '), end

% Take the ifft2 which completes the transform
fieldOut = ifft2(FT_Holo);
fieldOut = single(fieldOut*dx*dy);

% Crop due to padding
fieldOut = fieldOut(padding(1)+1:Nr-padding(1), padding(2)+1:Nc-padding(2));
if verbose, fprintf('Done.\n'), end