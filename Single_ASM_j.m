% Single hologram reconstruction at a single z or several positions using
% ASM. If the reconstruction distance is too high, a decimation is used
% until Fresnel reconstruction is allowed.
% Autofocusing algorithms : 2017 - Evaluation of refocus criteria for holographic particle imaging
% Distance of reconstruction from 2019 - Resolution and sampling analysis in digital in-line holography with spherical wave illumination
% z = (1/d - 1/L)^(-1) = Ld/(L-d)
% Made into function to use in JASMIN. Input a text file with information
% about the processing. See get_inputs help
% Effects of aliasing are considered. There is a maximum distance when
% using ASM (Fourier transfer) and minimum when using Fresnel convolution.
% Decimation by resampling and crop are allowed.
% Two windows give effective wavelength for reconstruction
% Always zero padding x2
%
% original_n_job_array : for using with script run_jobs_errored. This is
% the number of jobs in the original job array that produced nodes with
% errors. A list of jobids was created by join_Indices in the results
% folder with name not_found.mat.

function vout = Single_ASM_j(config_file,index,n_job_array,original_n_job_array)
if nargin == 1    % Single job case
   index = 1;
   n_job_array = 1;
else              % Job array case
   index = str2double(index);
   n_job_array = str2double(n_job_array);
end
if nargin == 4
   n_jobs_errored = n_job_array;
   n_job_array = str2double(original_n_job_array);
end
str_input = get_inputs(config_file);

image_path = str_input.image_path;
bgnd_path = str_input.bgnd_path;    % Empty to avoid
output_path = str_input.output_path;
zeval = str_input.zeval;
downsampling = str_input.downsampling;
crop_factor = str_input.crop_factor;
L = str_input.L;
make_movie_from_sweep = str_input.make_movie;
substract_mean_from_im = str_input.sub_mean_im;
num_pix_windows = str_input.size_windows;
auto_min_max_video = str_input.auto_min_max_video;
use_equivalent_wavelength = str_input.use_equiv_wavel;

% In case that we are dealing with jobs which node errored, transform
% job_index to the original job_index.
if nargin == 4
   loaded = load(fullfile(output_path,'not_found.mat'));
   if n_jobs_errored ~= length(loaded.not_found)
      error('The amount of detected failed jobs is different from the amount of jobs in the script array. Please correct to %d.', length(loaded.not_found))
   end
   index = loaded.not_found(index);
   fprintf('Redo job due to node fail.')
end
fprintf('Working with config file: %s. Job: %d/%d.\n',config_file,index,n_job_array)

movie_output_name = "movie z.avi";

% Name of file with movie output
min_val_movie = 20;
max_val_movie = 100;

dx0 = str_input.dx;
dy0 = str_input.dy;
wavelength = str_input.wavelength;

n1 = str_input.n1;
n2 = str_input.n2;
t1 = str_input.t1;
t2 = str_input.t2;
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
if ~isempty(str_input.z_autofoc_ini)
%     z_autofoc = linspace(str_input.z_autofoc_ini,str_input.z_autofoc_end,str_input.z_autofoc_n);
    [z_autofoc, step_to_resolution_ratio] = sample_z_sph(str_input.z_autofoc_ini,str_input.z_autofoc_end,str_input.z_autofoc_n,L,D,wavelength_equivalent);
    fprintf('The step-to-resolution ratio is %.2f.\n', step_to_resolution_ratio);
elseif str_input.z_resol > 0
    z_autofoc = (-str_input.z_span/2:str_input.z_resol:str_input.z_span/2)' + str_input.z_interest;
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
   
   if make_movie_from_sweep && n_job_array == 1
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
      fprintf('Propagation for index %d / %d : ',kz,Nz)
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
      % Movie maker
      if make_movie_from_sweep && n_job_array == 1
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
   end
   
   fprintf('Saving indices... ')
   clear('A', 'prevA', 'Psi', 'kernel_Fresnel', 'frame')
   all_indices = {Entropy, L1, Gradient, Laplacian, Variance, Tamura, CC, Gini}; % RC, 
   all_indices = cellfun(@(A) A./max(A),all_indices,'UniformOutput',false);
   if ~isfolder(output_path)
      mkdir(output_path);
   end
   if n_job_array == 1
      save(fullfile(output_path,'Indices.mat'),'all_indices', 'L', 'D', 'Nr0', 'Nc0', 'autof_Ws', 'z_autofoc')
   else
      save(sprintf(fullfile(output_path,'Indices_%04d.mat'),index),'all_indices', 'L', 'D', 'Nr0', 'Nc0', 'autof_Ws', 'z_autofoc')
   end
   fprintf('Done\n')
   
%     fig_indices = figure; plot(z_autofoc, all_indices)
%     legend('Entropy', 'L1', 'Gradient', 'Laplacian', 'Variance', 'Tamura', 'CC', 'RC', 'Gini')
%     saveas(fig_indices,'.\results\Indices.fig')
   
   if make_movie_from_sweep && n_job_array == 1
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

Nzi = length(z_recon);
for kz = 1:Nzi
   fprintf('kz = %d, index = %d, JOB_ID = %d. ',kz, ceil(kz/Nzi*n_job_array), index)
   if ceil(kz/Nzi*n_job_array) == index
      fprintf('Creating image %d. ', kz)
      down_required = abs(z_recon(kz)) / min(2*Nr0*dy0^2,2*Nc0*dx0^2) *wavelength_equivalent;   % Zero padding considered
      if down_required > downsampling
         warning('Occurring aliasing not being treated in z position %d of the output images.', kz)
      end
      fieldOut = easy_prop(k,z_recon(kz),root,FT_Holo,padding);
      fprintf('Propagation %d done.\n', kz)
      
      fieldOut = abs(fieldOut);
      fieldOut = uint8(fieldOut / max(fieldOut(:)) * 255.99);
      imwrite(fieldOut,fullfile(output_path,sprintf('im_output_%02d.png',kz)),'png')
      fprintf('Image %d saved.\n',kz)
      % Save im result in amplitude and phase
      % save('output_Single_ASM','fieldOut','Metrics','str_input')
   else
      fprintf('Not mine.\n')
   end
end
vout = 1;
fprintf('All done.\n')