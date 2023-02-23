function FT_Fresnel = get_Fresnel_FT(dx,dy,Nc,Nr,wavel,z)
k = 2*pi/wavel;

x = - Nc/2:Nc/2-1;  % xs and ys
y = - Nr/2:Nr/2-1;
[xx,yy] = meshgrid(x,y);        % root in phase multiplier
clear('x', 'y');

exp(1j*k*z)/(1j*wavel*z) * exp(1j*pi/wavel/z)(x^2+y^2)\right) \mathfrak{F}\left[U_o(x,y)\exp\left(\frac{i\pi}{\lambda z}(x^2+y^2))\right]