% Simulation of holograms obtained for a sample with a grid of dots as the
% one from Thorlabs R2L2S3P1.
% Using ASM. Planar wave.
% Axial resolution lambda/NA^2 : 2015 - Practical algorithms for simulation and reconstruction of digital in-line holograms
% Effects of aliasing are considered. There is a maximum distance when
% using ASM (Fourier transfer).
% Always zero padding x2

clear variables
close all
mm = 1e-3;
um = 1e-6;
nm = 1e-9;
LogicalStr = {'false', 'true'};

% Sample generation properties
dot_diam = 62.5*um;
dot_spacing = 125*um;
size_im = 10*mm;
dx = dot_diam/12;          % Size of pixel
N = 2*ceil(size_im/dx/2);  % Amount of pixels in a row/column
dy = dx;

output_path = './results/PeriodicSample';
log_file = fullfile(output_path, 'log.txt');
if ~isfolder(output_path)
   mkdir(output_path);
end

method = 'ASM';
wavelength = 355*nm;
% Distances of reconstruction
zmin = 1*mm;
zmax = 0.99 * 2*N*dx^2 / wavelength;   % Zero padding considered
Nz = 50;
zs = linspace(zmin,zmax, Nz);

substract_mean_from_im = false;

% Header of log file
logID = fopen(log_file,'w');
fprintf(logID, "Simulated dot pattern.\n\nDot diameter = %.2f um\nDot spacing = %.2f um\nResolution (dx) = %.2f um\nImage size = %dx%d\n\n",...
   dot_diam/um, dot_spacing/um,dx/um,N,N);
fprintf(logID, "Holographic properties\n\nWavelength = %.1f nm\nzs = %.2f mm : %.2f mm\nNumber of simulated distances = %d\nSubstraction of mean value before propagation : %s\n\n",...
   wavelength/nm, zmin/mm, zmax/mm, Nz, LogicalStr{substract_mean_from_im+1});
fclose(logID);

k = 2*pi / wavelength;

% Sample generation by oversampling
im = [];
im_file = fullfile(output_path,'im_simul.mat');
if exist(im_file,'file')
   s = load(im_file);
   if s.dot_diam == dot_diam && s.dot_spacing == dot_spacing && s.size_im == size_im && s.dx == dx
      im = s.im;
      logID = fopen(log_file); fprintf(logID, "Image loaded.\n"); fclose(logID);
      clear s
   end
end
if isempty(im)
   overs = 4;
   x = (1:overs*N)*dx/overs;
   dist_v2 = mod(x,dot_spacing)-dot_spacing/2;
   dist_v2 = dist_v2.^2;
   dist_im = sqrt(dist_v2'+dist_v2);
   im = dist_im > dot_diam/2;
   im = blockproc(im, [overs, overs], @(x) mean(x.data(:)));
   logID = fopen(log_file); fprintf(logID, "Image generation ready.\n"); fclose(logID);
   save(im_file,'im','dot_diam','dot_spacing','size_im','dx')
end

% Pad to double size
if substract_mean_from_im
   im = im - mean(im(:));
end
padding = [N/2, N/2];
if substract_mean_from_im
   im = padarray(im, padding);
else
   im = padarray(im, padding, 'replicate');
end
Nr = N*2; Nc = N*2;

detected_original_period = im_peak_period(im,dx);

% Prepare FFT2
FT_Holo = fftshift(fft2(im));
clear im

% Prepare kernel of ASM
root = get_ASM_root(dx,dy,Nc,Nr,wavelength);

logID = fopen(log_file); fprintf(logID,'Preparations for ASM ready.\n'); fclose(logID);

detected_periods = zeros(Nz,1);
wbar = waitbar(0, 'Processing z sweep...');
for kz = 1:Nz
   % Prepare this z
   z_for = zs(kz);
   % Check if propagation method and kernel are okay
   switch method
      case 'ASM'
         Psi = easy_prop(k, z_for, root, FT_Holo, padding, false); % Last input is Verbosity
      case 'Fresnel'
         [Psi, kernel_Fresnel] = propagate_Fresnel(dx,dy,wavelength_equivalent,z_for,FT_Holo,padding,kernel_Fresnel);
   end
   A = double(abs(Psi));
   clear Psi
   
   detected_periods(kz) = im_peak_period(A,dx);
   waitbar(kz/Nz,wbar);
end
close(wbar)

clear('A')

save(fullfile(output_path,'DetectedPeriods.mat'),'detected_periods')
logID = fopen(log_file);fprintf(logID, 'All done.\n'); fclose(logID);
plot(zs,detected_periods)