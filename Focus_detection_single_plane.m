% Single hologram reconstruction at several positions using ASM.
% If the reconstruction distance is too high, a decimation is used
% until Fresnel reconstruction is allowed.
% Autofocusing algorithms : 2017 - Evaluation of refocus criteria for holographic particle imaging
% Distance of reconstruction from 2019 - Resolution and sampling analysis in digital in-line holography with spherical wave illumination
% z = (1/d - 1/L)^(-1) = Ld/(L-d)
% Axial resolution lambda/NA^2 : 2015 - Practical algorithms for simulation and reconstruction of digital in-line holograms
% Effects of aliasing are considered. There is a maximum distance when
% using ASM (Fourier transfer) and minimum when using Fresnel convolution.
% There is an option of using band limited ASM.
% Decimation by resampling and crop are allowed.
% Two windows give effective wavelength for reconstruction
% Always zero padding x2
% Continue analysis with `analysis_Indices_focus_detection`

clear variables
close all
mm = 1e-3;
um = 1e-6;
nm = 1e-9;

images_path = 'D:\HoloICE\Holograms\Cold room tests\2022_05_11 Calibration\Grid';
format_images = 'tiff';
% Background image information
bgnd_path = []; % 'D:\HoloICE\Holograms\Pollen\Photron\ash_inv3Ms_maxpower_01\Calibration\bgnd\Camera_1_C001H001S0023_noHist.png'; % Background image full path or [].
get_ranges = true;   % Reduced sampling in z to know the z range in which to propagate. YOU NEED TO FIX zMin, zMax.

output_path = fullfile(images_path,'results\');
downsampling = 1; % It may produce aliasing effect. Jan suggested to use crop
crop_factor = 0.25;

zeval = []*mm;  % Set of visual inspection distances
zMin = [20*mm, 20*mm];   % Start and end values of focus detection range (based on order of images filenames.
zMax = [300*mm, 300*mm];
dx0 = 2.2 *um;  % Photron 20 um, Optronis 8 um, SVS Vistek 2.2 um 
dy0 = 2.2 *um;
wavelength = 355 *nm;
if get_ranges
   K_axial_sampling = 15;   % z sampling step-to-resolution ratio
else
   K_axial_sampling = 1;
end

% Autofocus parameter
num_pix_windows = [20, 100];

verbose = false;
substract_mean_from_im = false;
plot_detecting_focus = true;

% Get image names
fileList = dir(fullfile(images_path, sprintf('*.%s', format_images)));

Nim = length(fileList);
z_focus = zeros(Nim,1);
zMins = linspace(zMin(1), zMin(2), Nim);
zMaxs = linspace(zMax(1), zMax(2), Nim);
if get_ranges
   range_for = [1, Nim];
else
   range_for = 1:Nim;
end
for kim = range_for
   this_indices_file = fullfile(output_path,sprintf('Indices%02d.mat',kim));
   if exist(this_indices_file,'file')
      load(this_indices_file)
   else
      image_path = fullfile(images_path, fileList(kim).name);
      
      % Read images
      im = imread(image_path);
      if ~isempty(bgnd_path)
         im_bgnd = imread(bgnd_path);
      else
         im_bgnd = [];
      end
      % Using only good rows, or fixed image
%       [gr_start, gr_end, im1_corrected] = detect_good_row_start_end(im,im_bgnd);
      gr_start = 1; gr_end = size(im, 1); im1_corrected = [];
      if gr_start ~= 1 || gr_end ~= size(im, 1)
         fprintf('Using data from row %d to %d\n.', gr_start, gr_end);
      end
      if isempty(im1_corrected)
         im = im(gr_start:gr_end,:);
      else
         im = im1_corrected;
         clear im1_corrected
      end
      if ~isempty(im_bgnd)
         im_bgnd = im_bgnd(gr_start:gr_end,:);
         % Apply background
         im = single(im);
         im_bgnd = single(im_bgnd);
         im = (im - im_bgnd) ./ im_bgnd;
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
      
      k = 2*pi / wavelength;
      
      current_downsampling = downsampling;
      method = 'ASM';
      % Downsample
      [im, dx, dy, Nc, Nr] = downsample_image(im,downsampling,dx0,dy0);
      wanted_size = [Nr, Nc];
      
      % The z distances of reconstruction are based on a collimated illumination
      % beam. They are calculated as having constant step to axial resolution
      % ratio.
      zs = axial_sampling(zMins(kim), zMaxs(kim), [Nc, Nr], [dx dy], K_axial_sampling, wavelength);
      fprintf('Steps from %.2f um to %.2f um\n', (zs(2)-zs(1))*1e6, (zs(end)-zs(end-1))*1e6)
      
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
      meanim = mean(im(:));
      if substract_mean_from_im
         im = im - meanim;
      end
      im = im / meanim;
      padding = [Nr/2, Nc/2];
      if substract_mean_from_im
         im = padarray(im, padding);
      else
%          im = padarray(im, padding, 'replicate');
         im = padarray(im, padding, 1);
      end
      Nr = Nr * 2; Nc = Nc * 2;
      
      % Condition of suitability. Assuming object 20th of FOV
      l_object = max(Nc*dx, Nr*dy) /20;
      z_suit_Fresnel = 10*(pi * l_object^4 / 64 / wavelength)^(1/3);  % Minimum required distance
      % Condition of numeric implementation (dx=dy)
      z_imp_Fresnel = sqrt(Nr^2+Nc^2)*dx^2/wavelength;
      z_min_Fresnel = max(z_suit_Fresnel,z_imp_Fresnel);
      
      % Prepare FFT2
      FT_Holo = fftshift(fft2(im));
      clear im
      
      if verbose
         fprintf('FFT ready.\nResolution of images : %d x %d\nSize : %d MB\n',Nr,Nc,Nr*Nc*8/1024/1024); % Single precision (4 bytes) complex (x2)
      end
      
      % Prepare kernel of ASM
      root = get_ASM_root(dx,dy,Nc,Nr,wavelength);
      
      if verbose
         fprintf('root ASM ready.\n')
      end
      
      % Sweep z for autofocusing
      Nz = length(zs);
      kernel_Fresnel = [];
      
      % Initialize metrics
      StdAmpGrad = zeros(Nz, NW, 2);    % max
      StdComplexGrad = zeros(Nz, NW, 2);   % max
      
      EXP_MAX_MIN = 4;
      % Functions
      funs = autofocusing_funs();
      
      msg = 'Processing z sweep...';
      wbar = waitbar(0, msg);
      time_ini = clock;
      last_refresh = clock;
      refreshing_time = 2;    % seconds
      for kz = 1:Nz
         % Prepare this z
         z_for = zs(kz);
         % Check if propagation method and kernel are okay
         if strcmp(method,'ASM')
            down_required = abs(z_for) / min(2*Nr0*dy0^2,2*Nc0*dx0^2) *wavelength;   % Zero padding considered
            if down_required > current_downsampling || z_for > z_min_Fresnel
               fprintf('Modification of propagation required.')
               % Check if Fresnel propagation is available
               if z_for > z_min_Fresnel
                  method = 'Fresnel';
                  current_downsampling = downsampling;
                  warning('Propagation method was changed to Fresnel convolution in reconstruction %d.',kz);
               else
                  current_downsampling = ceil(down_required);
                  warning('Downsampling factor was changed to avoid aliasing in reconstruction %d to value = %d.',kz,current_downsampling);
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
                     root = get_ASM_root(dx,dy,Nc,Nr,wavelength);
                  end
                  fprintf(' new root ready.')
               end
               fprintf('Modification of propagation done.\n')
            end
         end
         switch method
            case 'ASM'
               Psi = easy_prop(k, z_for, root, FT_Holo, padding, verbose);
            case 'Fresnel'
               [Psi, kernel_Fresnel] = propagate_Fresnel(dx,dy,wavelength,z_for,FT_Holo,padding,kernel_Fresnel);
         end
         % Upsample if necessary
         if current_downsampling > downsampling
            Psi = upsample_result(Psi,current_downsampling/downsampling,wanted_size);
         end
         A = double(abs(Psi));
         % clear Psi % Comment if StdComplexGrad or RC is used
         
         if verbose, fprintf('Autofocusing : '), end
         for kW = 1:NW
            autof_W = autof_Ws(kW);
            % Std Amplitude
            B = blockproc(A,[autof_W autof_W],funs.StdGrad,'DisplayWaitbar',false);
            StdAmpGrad(kz, kW, 1) = max(B(:));
            StdAmpGrad(kz, kW, 2) = sum(B(:).^EXP_MAX_MIN);
            if verbose, fprintf('.'), end
            % Std Complex
            B = blockproc(Psi,[autof_W autof_W],funs.StdGrad,'DisplayWaitbar',false);
            StdComplexGrad(kz, kW, 1) = max(B(:));
            StdComplexGrad(kz, kW, 2) = sum(B(:).^EXP_MAX_MIN);
            if verbose, fprintf('.|'), end
         end
         if verbose, fprintf(' Indices found.\n'), end
         
         if etime(clock, last_refresh) > refreshing_time
            last_refresh = clock;
            elapsed_time = etime(clock, time_ini);
            remaining_time = (Nz - kz) * elapsed_time / kz;
            if remaining_time > 120
               msg = sprintf('Processing z sweep... %d minutes remaining.', round(remaining_time/60));
            else
               msg = sprintf('Processing z sweep... %d seconds remaining.', round(remaining_time));
            end
         end
         waitbar(kz/Nz,wbar, msg);
      end
      close(wbar)
      
      fprintf('Saving indices... ')
      clear('A', 'Psi', 'kernel_Fresnel', 'frame')
      all_indices = {StdAmpGrad, StdComplexGrad};
      all_indices = cellfun(@(A) A./max(A),all_indices,'UniformOutput',false);
      legend = {'StdAmpGrad', 'StdComplexGrad'};
      
      if ~exist(output_path, 'dir')
         mkdir(output_path);
      end
      save(this_indices_file,'all_indices', 'legend', 'autof_Ws', 'zs')
      %    fprintf('Done\n')
      %    plot_indices([output_path 'Indices.mat'],zs,autof_Ws,'sum',.5,legend)   % Sums results over all windows
   end
   
   n_win = size(all_indices{1},2);
   select_sizes = n_win-1:n_win;   % Don't use small windows
   z_focus(kim) = detect_focus_z(all_indices,zs,autof_Ws,select_sizes,[],plot_detecting_focus,output_path,num2str(kim));
   fprintf('Image %d of %d processed.\n', kim, Nim)
end
if ~get_ranges
   save(fullfile(output_path,'z_positions.mat'), 'z_focus', 'fileList')
end