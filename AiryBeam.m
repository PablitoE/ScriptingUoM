% Simulation of Airy beams to check on the influence of resolution and
% pixel size
% 2021 - Optical trapping with non-diffracting Airy beams array using a
% holographic optical tweezers
clear variables
close all
mm = 1e-3;
um = 1e-6;
nm = 1e-9;

L = 15.4*mm;      % Similar to Thorlabs HD2 4k, Holoeye GAEA, Hamamatsu
dp = 8*um; % Thorlabs HD2, 6.4 HD1, 3.74 4k, Hamamatsu 12.5, Holoeye Pluto 8, LUNA 4.5, GAEA 3.74, Meadowlark 9.2

Np = floor(L/dp);    % Number of pixels
Nz = 1000;           % Number of z positions to analyse
Lz = 1000*mm;        % Distance of propagation in z to analyse
wavel = 532*nm;      % Wavelength
x0 = .1*mm;          % Scale factor

nu = 2;            % Initial orientation
a = 0.07;            % Convergence constant

do_1D = false;
do_beams = true;

%%%%%%%%%%%%%%%%% 1D analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%
AiryBeamFun = @(ms,mz) airy(ms - mz.^2/4 - nu*mz + 1i*a*mz).* exp(a*(ms-mz.^2/2-nu*mz)) .*...
   exp(1i*(-mz.^3/12 + (a^2-nu^2+ms).*mz/2 + nu*ms - nu*mz.^2/2));

if do_1D
   s = (linspace(-L/2,L/2,Np)+L/4)/x0;
   z = linspace(0,Lz/2/pi*wavel/x0^2,Nz);
   
   [ss,zz] = meshgrid(s,z);
   Psi = AiryBeamFun(ss,zz);
   Psi_a = abs(Psi);
   imagesc(Psi_a)
end
%%%%%%%%%%%%%%%%%%% Airy beams array %%%%%%%%%%%%%%%%%%%%%
if do_beams
   Nab = 4;          % Number of Airy beams
   xp = 1*mm; % L/4;         % Displacement of each beam from the center pixel
   thetas = 2*(0:Nab-1)*pi/Nab;  % Angles of beams in a circle with radius xp
   focal_point = 4*pi/wavel*x0^2*(sqrt(nu^2+xp/x0)-nu);  % Added ^2 for units
   max_z_sim = 1.5*focal_point;
   % Prepare some transversal images and one longitudinal image (integrated)
   Nim = 10;

   thetas_ = reshape(thetas,[1,1,Nab]);
   dxs = xp * cos(thetas_);
   dys = xp * sin(thetas_);
   x = linspace(-L/2,L/2,Np);
   [xx,yy] = meshgrid(x);
   sjx = (xx.*cos(thetas_)-yy.*sin(thetas_)+xp)/x0;
   sjy = (xx.*sin(thetas_)+yy.*cos(thetas_)+xp)/x0;
   % Ims z
   true_zim = linspace(0,max_z_sim,Nim);
   zim = linspace(0,max_z_sim/2/pi*wavel/x0^2,Nim);
   Psi_imz = zeros(Np,Np,Nim);
   for jz = 1:Nim
      for ja = 1:Nab
         Psix = AiryBeamFun(sjx(:,:,ja),zim(jz));
         Psiy = AiryBeamFun(sjy(:,:,ja),zim(jz));
         Psi_imz(:,:,jz) = Psi_imz(:,:,jz) + Psix.*Psiy;
      end
   end
   
   for jz =1:Nim
      imagesc(abs(Psi_imz(:,:,jz)))
      title(sprintf('z = %.2f mm (focus at %.2f mm.',true_zim/mm,focal_point/mm))
      pause
   end
end