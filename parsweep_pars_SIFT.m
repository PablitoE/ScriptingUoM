% Analysis of Indices.mat, results from Focus_detection_single_plane.
clear variables
close all
path_to_file = 'D:\HoloICE\Scripting\results\PhotronUSAFCali\spline_max\z_positions.mat';
% If given this folder, the images are extracted from there
ims_folder = 'D:\HoloICE\Holograms\Glass beads\Photron - 2um\20210812\Calibration USAF';

load(path_to_file)

% plot_indices(path_to_file,zs,autof_Ws,'sum',[],legend)   % Filtered
if exist('results\sweep_pars_SIFT_out.mat','file')
   load('results\sweep_pars_SIFT_out.mat')
   Npars = Np * Ne * Nt * Nh;
else
   sector_bb = [257, 768, 256, 768];
   dx = 20e-6;
   wavelength = 532e-9;
   
   Np = 20; Ne = 15; Nt = 10; Nh = 5;
   Npars = Np * Ne * Nt * Nh;
   peaks_th = logspace(log10(3), log10(50), Np);
   edges_th = logspace(log10(3), log10(30), Ne);
   thrs_ubc = logspace(log10(0.9), log10(4), Nt);
   hp_freqs = logspace(log10(0.02), log10(0.2), Nh);
   
   MAX_TOL_ERR = 3; % Maximum tolerated error when detecting outliers in matched features by comparing with expected transformed location of initial features.
   
   Nim = length(fileList);
   
   mean_goods = zeros(Npars,Nim-1);
   std_mean_goods = zeros(Npars,Nim-1);
   max_goods = zeros(Npars,Nim-1);
   
   msg = 'Processing sweep...';
   % wbar = waitbar(0, msg);
   time_ini = clock;
   last_refresh = clock;
   refreshing_time = 5;    % seconds
   WaitMessage = parfor_wait(Npars, 'Waitbar', true);
   
   parfor ki = 1:Npars
      [kh, kt, ke, kp] = ind2sub([Nh,Nt,Ne,Np], ki);
      peak_thresh = peaks_th(kp);
      edge_thresh = edges_th(ke);
      thresh_ubc = thrs_ubc(kt);
      hp_freq = hp_freqs(kh);
      
      % Create highpass filter
      [f1,f2] = freqspace(21,'meshgrid');
      r = sqrt(f1.^2 + f2.^2);
      Hd = ones(21);
      Hd(r<hp_freq) = 0;
      win = fspecial('gaussian',21,2);
      win = win ./ max(win(:));
      h = fwind2(Hd,win);
      % clear f1 f2 r Hd win
      
      slice_mean_goods = zeros(Nim-1,1);
      slice_std_mean_goods = zeros(Nim-1,1);
      slice_max_goods = zeros(Nim-1,1);
      
      % Deal with first image
      if exist(ims_folder,'dir')
         image_path = fullfile(ims_folder, fileList(1).name);
      else
         image_path = fullfile(fileList(1).folder, fileList(1).name);
      end
      % Read images
      im = imread(image_path);
      Psi = hologram_reconstruction_at_z(im, z_focus(1), wavelength, dx);
      aPsi_prev = abs(Psi);
      aPsi_prev = filter2(h, aPsi_prev);
      aPsi_min = min(aPsi_prev(:));
      aPsi_max = max(aPsi_prev(:));
      aPsi_prev = single((aPsi_prev-aPsi_min)/(aPsi_max-aPsi_min)*255);
      [f_prev,d_prev] = vl_sift(aPsi_prev, 'PeakThresh', peak_thresh, 'edgethresh', edge_thresh);
      
      for kim = 2:Nim
         if exist(ims_folder,'dir')
            image_path = fullfile(ims_folder, fileList(kim).name);
         else
            image_path = fullfile(fileList(kim).folder, fileList(kim).name);
         end
         % Read images
         im = imread(image_path);
         Psi = hologram_reconstruction_at_z(im, z_focus(kim), wavelength, dx);
         aPsi = abs(Psi);
         % clear Psi
         aPsi = filter2(h, aPsi);
         
         % SIFT vlFeat
         aPsi = single((aPsi-aPsi_min)/(aPsi_max-aPsi_min)*255);
         [f,d] = vl_sift(aPsi, 'PeakThresh', peak_thresh, 'edgethresh', edge_thresh);
         [matches, scores] = vl_ubcmatch(d_prev, d, thresh_ubc);
         matchedDistortedXY = [f(1,matches(2,:)); f(2,matches(2,:))]';
         matchedOriginalXY = [f_prev(1,matches(1,:)); f_prev(2,matches(1,:))]';
         [tformBackToOriginal,inlierIdx,err_struct] = my_estimateGeometricTransform2D(matchedDistortedXY,...
            matchedOriginalXY,'similarity',MAX_TOL_ERR,false);
         % show_match_sift(aPsi_prev,aPsi,matches,f_prev,f,inlierIdx,tformBackToOriginal);
         f_prev = f;
         d_prev = d;
         
         slice_mean_goods(kim-1) = err_struct.mean_good;
         slice_std_mean_goods(kim-1) = err_struct.std_good / sqrt(err_struct.n_good);
         slice_max_goods(kim-1) = err_struct.max_good;
         aPsi_prev = aPsi;
      end
      mean_goods(ki,:) = slice_mean_goods;
      std_mean_goods(ki,:) = slice_std_mean_goods;
      max_goods(ki,:) = slice_max_goods;
      
      WaitMessage.Send;
   end
   % close(wbar)
   WaitMessage.Destroy
   save('results\sweep_pars_SIFT_out.mat','mean_goods','std_mean_goods','max_goods','peaks_th','edges_th','thrs_ubc','hp_freqs','path_to_file','Nh','Nt','Ne','Np')
end
max_mean = max(mean_goods,[],2);
max_std_mean = max(std_mean_goods,[],2);
max_max = max(max_goods,[],2);

fprintf('Product of Mean, Std_Mean and Max:\n')
prod_maxs = max_mean .* max_std_mean .* max_max;
[~, i_prod] = min(prod_maxs);
[kh, kt, ke, kp] = ind2sub([Nh,Nt,Ne,Np], i_prod);
peak_thresh = peaks_th(kp)   % 15
edge_thresh = edges_th(ke)   % 3
thresh_ubc = thrs_ubc(kt)    % 2
hp_freq = hp_freqs(kh)       % 0.02

fprintf('Product of positions of Mean, Std_Mean and Max:\n')
is_m = {[],[],[]};
positions = is_m;
prod_positions = zeros(Npars,1);
[~,is_m{1}] = sort(max_mean);
[~,is_m{2}] = sort(max_std_mean);
[~,is_m{3}] = sort(max_max);
for kis = 1:3
   aux = [is_m{kis}, (1:Npars)'];
   aux = sortrows(aux);
   positions{kis} = aux(:,2);
   prod_positions = positions{kis} .* prod_positions;
end
[~, i_prod_positions] = min(prod_positions);
[kh, kt, ke, kp] = ind2sub([Nh,Nt,Ne,Np], i_prod_positions);
peak_thresh = peaks_th(kp)   % 3
edge_thresh = edges_th(ke)   % 3
thresh_ubc = thrs_ubc(kt)    % 0.9
hp_freq = hp_freqs(kh)       % 0.02

fprintf('Product of Mean and Max:\n')
prod_maxs = max_mean .* max_max;
[~, i_prod] = min(prod_maxs);
[kh, kt, ke, kp] = ind2sub([Nh,Nt,Ne,Np], i_prod);
peak_thresh = peaks_th(kp)   % 15
edge_thresh = edges_th(ke)   % 3
thresh_ubc = thrs_ubc(kt)    % 2
hp_freq = hp_freqs(kh)       % 0.02

fprintf('Product of positions of Mean, Std_Mean and Max:\n')
prod_positions = positions{1} .* positions{3};
[~, i_prod_positions] = min(prod_positions);
[kh, kt, ke, kp] = ind2sub([Nh,Nt,Ne,Np], i_prod_positions);
peak_thresh = peaks_th(kp)   % 15
edge_thresh = edges_th(ke)   % 3
thresh_ubc = thrs_ubc(kt)    % 2
hp_freq = hp_freqs(kh)       % 0.02