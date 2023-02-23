% Single hologram reconstruction at a single z or several positions using
% ASM. If the reconstruction distance is too high, a decimation is used
% until Fresnel reconstruction is allowed.
% Autofocusing algorithms : 2017 - Evaluation of refocus criteria for holographic particle imaging
% Distance of reconstruction from 2019 - Resolution and sampling analysis in digital in-line holography with spherical wave illumination
% z = (1/d - 1/L)^(-1) = Ld/(L-d)
% Axial resolution lambda/NA^2 : 2015 - Practical algorithms for simulation and reconstruction of digital in-line holograms
% Effects of aliasing are considered. There is a maximum distance when
% using ASM (Fourier transfer) and minimum when using Fresnel convolution.
% Decimation by resampling and crop are allowed.
% Two windows give effective wavelength for reconstruction
% Always zero padding x2

clear variables
close all
mm = 1e-3;
um = 1e-6;
nm = 1e-9;

image_path = 'C:\Users\j05625pe\Documents\Work\Holograms\Spherical R2L2S3P1\20201117\Background00_4pulses_down4.png';   
% Background image information
bgnd_path = []; % 'C:\Users\j05625pe\Documents\Work\Holograms\Spherical R2L2S3P1\20201030\Background00.tif'; % Background image. Good rows from 73

output_path = './results/';
zeval = [10.92, 18.98, 21.89, 28.19, 37.93]*mm; 
% [20.4 29.4 42.7 50.8 55 57.3 63.7 64.7 66.8 67.5]*mm;    % [24.65, 29.50, 42.1, 69]*mm; % 20201117\Background00_4pulses_down4.png, L = 80 mm , JASMIN and HoloViewer
% [26.06 26.76 26.96 29.42]*mm; % (gradient, laplacian, variance, RC) : 20201026\Im06
% [19.78, 27.56, 41.28, 43.43, 48.3, 50.23, 54.23]*mm; % % 20201026\Im00.tif L_equiv % Output images
% [21.77, 30.23, 44.97, 47.05, 52.43, 54.47] *mm; % 20201026\Im00.tif wavel_equiv % Output images
downsampling = 1; % It may produce aliasing effect. Jan suggested to use crop
crop_factor = 1;
L = 90 *mm; % Distance from sensor to focus point. Not 60.1 mm
% Distances of reconstruction
z_autofoc = []; % linspace(8*mm,21*mm, 50); % 
z_resol = 0; % 30*um;   % En 0 anula el sweep
z_span = 1*mm;
z_interest = [11, 19]*mm; % [22, 28, 32.6, 37.8, 45.7]*mm;
% [21.3, 29.5, 55, 58.1, 65, 67.9]*mm; % 20201117\Background00 L aventurado 80mm
% 48.3*mm; % 20201026\Im03.tif dist equiv
% z_interest = [21.9, 30.2, 45, 47.1, 52.4, 54.4] * mm; % 20201026\Im00.tif wavel_equiv
% z_interest = [19.75, 27.5, 41.35, 43.5, 48.27, 50.22, 54.3] * mm; % 20201026\Im00.tif dist equiv
% Autofocus parameter
num_pix_windows = [7, 100];

save_box_question = false; % Prompt a question to get x and y ranges
make_movie_from_sweep = false;
substract_mean_from_im = true;
save_output_images = true;    % Not shown if saved (only amplitude)
use_equivalent_wavelength = false; % If false, the total L will be modified to use real wavelength (n=1)

% Name of file with movie output
auto_min_max_video = true; % Manual values will be overwritten if true
min_val_movie = 0;
max_val_movie = 100;
movie_output_name = '.\results\movie z.avi';

dx0 = 2.2 *um;
dy0 = 2.2 *um;
wavelength = 355 *nm;

% Fused silica refractive index
% nFS = sqrt(1+.6961663 * wavelength^2/(wavelength^2-.0684043^2) + .4079426 * ...
%    wavelength^2/(wavelength^2-.1162414^2) + .8974794 * wavelength^2/...
%    (wavelength^2-9.896161^2));
nFS = 1.4761; % Fused Silica at 355 nm
n1 = nFS;
t1 = 5 *mm;
n2 = nFS;
t2 = 12 *mm;

% Dealing with windows effects
if use_equivalent_wavelength
   wavelength_equivalent = wavelength * L / (n1*t1+n2*t2 + (L-t1-t2));
else
   t1_eq = t1 / n1;
   t2_eq = t2 / n2;
   L = L - t1 + t1_eq - t2 + t2_eq;
   wavelength_equivalent = wavelength;
end

% Read images
im = imread(image_path);
if ~isempty(bgnd_path)
   im_bgnd = imread(bgnd_path);
else
   im_bgnd = [];
end
% Using only good rows, or fixed image
[gr_start, gr_end, im1_corrected] = detect_good_row_start_end(im,im_bgnd);
if isempty(im1_corrected)
   im = im(gr_start:gr_end,:);
else
   im = im1_corrected;
   clear im1_corrected
end
if ~isempty(im_bgnd)
   im_bgnd = im_bgnd(gr_start:gr_end,:);
   % Apply background
   im = im ./ im_bgnd;
end
clear im_bgnd

[Nr, Nc] = size(im); % [Ny, Nx]
% Cropping
Nr_c = floor(Nr * crop_factor); Nc_c = floor(Nc * crop_factor); 
ir_c = floor(Nr/2- Nr_c/2)+1; ir_c = ir_c:min(ir_c+Nr_c,Nr);
ic_c = floor(Nc/2- Nc_c/2)+1; ic_c = ic_c:min(ic_c+Nc_c,Nc);
im = im(ir_c,ic_c);
im0 = uint8(im);  % Save image for future required downsamplings
[Nr0, Nc0] = size(im); % [Ny, Nx]
D = max(Nr0 * dy0, Nc0 * dx0);

% Propagation distances
z_autofoc = [];
if ~isempty(z_autofoc_ini)
    [z_autofoc, step_to_resolution_ratio] = sample_z_sph(z_autofoc_ini,z_autofoc_end,z_autofoc_n,L,D,wavelength_equivalent);
    fprintf('The step-to-resolution ratio is %.2f.\n', step_to_resolution_ratio);
elseif z_resol > 0
    z_autofoc = (-z_span/2:z_resol:z_span/2)' + z_interest;
    z_autofoc = z_autofoc(:);
end
if any(z_autofoc >= L)
   warning("Some z positions from sweep were extracted from the analysis as being greater or equal to the equivalent distance of %.2f mm",L*1e3)
   z_autofoc = z_autofoc(z_autofoc < L);
end

% Reconstruction distances using plane wave
z_recon = zeval * L ./ (L-zeval);
z_recon_autofoc = [];
if ~isempty(z_autofoc)
   z_recon_autofoc = z_autofoc * L ./ (L-z_autofoc);
end

k = 2*pi / wavelength_equivalent;

% Required downsampling to avoid aliasing in the first reconstruction
% if isempty(z_autofoc)
%    z_check = z_recon(1);
% else
%    z_check = z_recon_autofoc(1);
% end
% down_required = abs(z_check) / min(2*Nr0*dy0^2,2*Nc0*dx0^2) *wavelength_equivalent;   % Zero padding considered
% if down_required > downsampling
%    downsampling = ceil(down_required);
%    warning('Downsampling factor was change to avoid aliasing in reconstruction to value = %.2f.',downsampling);
% end
current_downsampling = downsampling;
method = 'ASM';
% Downsample
[im, dx, dy, Nc, Nr] = downsample_image(im,downsampling,dx0,dy0);
wanted_size = [Nr, Nc];
% Prepare indices to extract a piece of the reconstruction selected by user
sbx = 0;
sby = 0;
if save_box_question
   imagesc(im)
   axis image
   [x, y] = ginput(2);
   close()
   x = floor(x); y = floor(y);
   sbx = min(x):max(x);
   sby = min(y):max(y);
end

% Get sizes of window analysis of autofocusing
if length(num_pix_windows) == 1
   autof_Ws = num_pix_windows;
elseif length(num_pix_windows) == 2 % Two inputs of the min and max amounts of windows in the smaller dimension (doubling size of window)
   min_autof_W = floor(min(num_pix_windows));
   max_autof_W = max(num_pix_windows);      % Possibly unreachable as integer
   n = floor(log(max_autof_W/min_autof_W)/log(2));
   n = 2.^(0:n);
   autof_Ws = min_autof_W * n;
   str_autof_Ws = sprintf('%d ',autof_Ws);
   fprintf('The window sizes are %s.\n',str_autof_Ws)
else
   error('Wrong input for window size for autofocusing.')
end
NW = length(autof_Ws);

% Pad to double size
if substract_mean_from_im
   im = im - mean(im(:));
end
padding = [Nr/2, Nc/2];
if substract_mean_from_im
   im = padarray(im, padding);
else
   im = padarray(im, padding, 'replicate');
end
Nr = Nr * 2; Nc = Nc * 2;

% Condition of suitability. Assuming object 20th of FOV
l_object = max(Nc*dx, Nr*dy) /20;
z_suit_Fresnel = 10*(pi * l_object^4 / 64 / wavelength_equivalent)^(1/3);  % Minimum required distance
% Condition of numeric implementation (dx=dy)
z_imp_Fresnel = sqrt(Nr^2+Nc^2)*dx^2/wavelength_equivalent;
z_min_Fresnel = max(z_suit_Fresnel,z_imp_Fresnel);

% Prepare FFT2
FT_Holo = fftshift(fft2(im));
clear im

fprintf('FFT ready.\nResolution of images : %d x %d\nSize : %d MB\n',Nr,Nc,Nr*Nc*8/1024/1024); % Single precision (4 bytes) complex (x2)

% Prepare kernel of ASM
root = get_ASM_root(dx,dy,Nc,Nr,wavelength_equivalent);

fprintf('root ASM ready.\n')

% Sweep z for autofocusing
if ~isempty(z_recon_autofoc)
   Nz = length(z_recon_autofoc);
   % Workable indices for this job in the array
   w_ind = 1:Nz;
   w_ind = w_ind(ceil(w_ind/Nz*n_job_array) == index);
   Nz = length(w_ind);
   kernel_Fresnel = [];

   if save_box_question
      saved_box = zeros(length(sby),length(sbx),Nz);
   end
   if make_movie_from_sweep
       writerObj = VideoWriter(movie_output_name);
       writerObj.FrameRate = 20;
       open(writerObj);
   end
   % Initialize metrics
   Entropy = zeros(Nz, NW, 2);    % max
   L1 = zeros(Nz, NW, 2);         % min
   Gradient = zeros(Nz, NW, 2);   % max
   Laplacian = zeros(Nz, NW, 2);  % max
   Variance = zeros(Nz, NW, 2);   % max
   Tamura = zeros(Nz, NW, 2);     % max
   CC = zeros(Nz, NW, 2);         % min
%    RC = zeros(Nz, NW, 2);         % max
   Gini = zeros(Nz, NW, 2);       % min
   EXP_MAX_MIN = 4;
   % Functions
   funs = autofocusing_funs();
   
   wbar = waitbar(0, 'Processing z sweep...');
   for kz = 1:Nz
      % Prepare this z
      z_for = z_recon_autofoc(w_ind(kz));
      % Check if propagation method and kernel are okay
      if strcmp(method,'ASM')
         down_required = abs(z_for) / min(2*Nr0*dy0^2,2*Nc0*dx0^2) *wavelength_equivalent;   % Zero padding considered
         if down_required > current_downsampling || z_for > z_min_Fresnel
            fprintf('Modification of propagation required.')
            % Check if Fresnel propagation is available
            if z_for > z_min_Fresnel
               method = 'Fresnel';
               current_downsampling = downsampling;
               warning('Propagation method was changed to Fresnel convolution in reconstruction %d.',w_ind(kz));
            else
               current_downsampling = ceil(down_required);
               warning('Downsampling factor was changed to avoid aliasing in reconstruction %d to value = %d.',w_ind(kz),current_downsampling);
            end
            new_dx = dx0*current_downsampling;  % Check the need of recalculating
            if dx ~= new_dx
               fprintf('Preparing downsampling.')
               [im, dx, dy, Nc, Nr] = downsample_image(im0,current_downsampling,dx0,dy0);
               % Pad to double size
               if substract_mean_from_im
                  im = im - mean(im(:));
               end
               padding = [Nr/2, Nc/2];
               im = padarray(im, padding); % , 'replicate');
               Nr = Nr * 2; Nc = Nc * 2;
               % Prepare FFT2
               FT_Holo = fftshift(fft2(im));
               clear im
               fprintf(' FFT ready.')
               % Prepare kernel of propagation
               if strcmp(method,'ASM')
                  root = get_ASM_root(dx,dy,Nc,Nr,wavelength_equivalent);
               end
               fprintf(' new root ready.')
            end
            fprintf('Modification of propagation done.\n')
         end
      end
      switch method
         case 'ASM'
            Psi = easy_prop(k, z_for, root, FT_Holo, padding);
         case 'Fresnel'
            [Psi, kernel_Fresnel] = propagate_Fresnel(dx,dy,wavelength_equivalent,z_for,FT_Holo,padding,kernel_Fresnel);
      end
      fprintf('Done.\n')
      % Upsample if necessary
      if current_downsampling > downsampling
         Psi = upsample_result(Psi,current_downsampling/downsampling,wanted_size);
      end
      A = double(abs(Psi));
      clear Psi % Comment if RC is used
      % Saved piece
      if save_box_question
         saved_box(:,:,kz) = A(sby,sbx);
      end
      % Movie maker
      if make_movie_from_sweep
         if kz == 1 && auto_min_max_video
            min_val_movie = max(0, min(A(:))-10);
            max_val_movie = min(255, max(A(:))*1.2);
         end
         frame = (A-min_val_movie)/(max_val_movie-min_val_movie);
         frame(frame>1) = 1;
         frame(frame<0) = 0;
         frame = im2frame(repmat(frame,[1 1 3]));
         writeVideo(writerObj, frame);
      end

      fprintf('Autofocusing : ')
      for kW = 1:NW
         autof_W = autof_Ws(kW);
         % Entropy (amplitude) : 54.8 mm 56 mm
         B = blockproc(A,[autof_W autof_W],funs.Entropy,'DisplayWaitbar',false);
         Entropy(kz, kW, 1) = max(B(:));
         Entropy(kz, kW, 2) = sum(B(:).^EXP_MAX_MIN);
         fprintf('.')
         % L1 norm
         B = blockproc(A,[autof_W autof_W],funs.L1,'DisplayWaitbar',false);
         L1(kz, kW, 1) = min(B(:));
         L1(kz, kW, 2) = sum(B(:).^(-EXP_MAX_MIN));
         fprintf('.')
         % Gradient
         B = blockproc(A,[autof_W autof_W],funs.Gradient,'DisplayWaitbar',false);
         Gradient(kz, kW, 1) = max(B(:));
         Gradient(kz, kW, 2) = sum(B(:).^EXP_MAX_MIN);
         fprintf('.')
         % Laplacian
         B = blockproc(A,[autof_W autof_W],funs.Laplacian,'DisplayWaitbar',false);
         Laplacian(kz, kW, 1) = max(B(:));
         Laplacian(kz, kW, 2) = sum(B(:).^EXP_MAX_MIN);
         fprintf('.')
         % Variance
         B = blockproc(A,[autof_W autof_W],funs.Variance,'DisplayWaitbar',false);
         Variance(kz, kW, 1) = max(B(:));
         Variance(kz, kW, 2) = sum(B(:).^EXP_MAX_MIN);
         fprintf('.')
         % Tamura
         B = blockproc(A,[autof_W autof_W],funs.Tamura,'DisplayWaitbar',false);
         Tamura(kz, kW, 1) =  max(B(:));
         Tamura(kz, kW, 2) =  sum(B(:).^EXP_MAX_MIN);
         fprintf('.')
         % Correlation coefficient
         if kz == 1
            CC(kz, kW, 1) = 0;
            CC(kz, kW, 2) = 0;
         else
            B = funs.CC(A, prevA, autof_W);
            CC(kz, kW, 1) = min(B(:));
            CC(kz, kW, 2) = sum(B(:).^(-EXP_MAX_MIN));
         end
         fprintf('.')
%          % RC method
%          B = blockproc(Psi,[autof_W autof_W],funs.RC,'DisplayWaitbar',false);
%          RC(kz, kW, 1) = max(B(:));
%          RC(kz, kW, 2) = sum(B(:).^EXP_MAX_MIN);
%          fprintf('.')
         % Gini
         B = blockproc(A,[autof_W autof_W],funs.Gini,'DisplayWaitbar',false);
         Gini(kz, kW, 1) = min(B(:));
         Gini(kz, kW, 2) = sum(B(:).^(-EXP_MAX_MIN));
         fprintf('|')
      end
      prevA = single(A);
      fprintf(' Indices found.\n')
      waitbar(kz/Nz,wbar);
   end
   close(wbar)

   fprintf('Saving indices... ')
   clear('A', 'prevA', 'Psi', 'kernel_Fresnel', 'frame')
   all_indices = {Entropy, L1, Gradient, Laplacian, Variance, Tamura, CC, Gini}; % RC,
   all_indices = cellfun(@(A) A./max(A),all_indices,'UniformOutput',false);
   % fig_indices = figure; plot(z_autofoc, all_indices)
   % legend('Entropy', 'L1', 'Gradient', 'Laplacian', 'Variance', 'Tamura', 'CC', 'RC', 'Gini')
%    saveas(fig_indices,'.\results\Indices.fig')

   if ~isfolder(output_path)
      mkdir(output_path);
   end
   save(fullfile(output_path,'Indices.mat'),'all_indices')
   fprintf('Done\n')
   if z_resol == 0
      plot_indices([output_path 'Indices.mat'],z_autofoc,autof_Ws,'sum',.5)   % Sums results over all windows
   else
      plot_indices([output_path 'Indices.mat'],z_autofoc,autof_Ws,'sum')   % Sums results over all windows
   end
   
   if make_movie_from_sweep
      close(writerObj);
   end
end

fprintf('Starting with particular images.\n')
% Propagation to desired z
if ~isempty(z_recon_autofoc)
   [im, dx, dy, Nc, Nr] = downsample_image(im0,downsampling,dx0,dy0);
   % Pad to double size
   if substract_mean_from_im
      im = im - mean(im(:));
   end
   padding = [Nr/2, Nc/2];
   im = padarray(im, padding); % , 'replicate');
   Nr = Nr * 2; Nc = Nc * 2;
   % Prepare FFT2
   FT_Holo = fftshift(fft2(im));
   clear im
   % Prepare kernel of propagation
   root = get_ASM_root(dx,dy,Nc,Nr,wavelength_equivalent);
end

wbar = waitbar(0, 'Processing output images...');
Nzi = length(z_recon);
for kz = 1:Nzi
   fprintf('Creating image %d at distance %.2f mm. ', kz, z_recon(kz)/mm)
   down_required = abs(z_recon(kz)) / min(2*Nr0*dy0^2,2*Nc0*dx0^2) *wavelength_equivalent;   % Zero padding considered
   if down_required > downsampling
      warning('Occurring aliasing not being treated in z position %d of the output images.', kz)
   end
   fieldOut = easy_prop(k,z_recon(kz),root,FT_Holo,padding);
   fprintf('Propagation %d done.\n', kz)
   if save_output_images
      fieldOut = abs(fieldOut);
      fieldOut = uint8(fieldOut / max(fieldOut(:)) * 255.99);
      imwrite(fieldOut,fullfile(output_path,sprintf('im_output_%02d.png',kz)),'png')
      fprintf('Image %d saved.\n',kz)
   else
      % Show result in amplitude and phase
      figure(kz), subplot(1,2,1)
      imagesc(abs(fieldOut))
      colormap('gray'), axis image
      subplot(1,2,2), imagesc(angle(fieldOut))
      axis image
   end
   waitbar(kz/Nzi,wbar);
end
close(wbar)
fprintf('All done.\n')