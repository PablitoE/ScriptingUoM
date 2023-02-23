% Analysis of Indices.mat, results from Focus_detection_single_plane.
clear variables
close all
path_to_file = 'D:\HoloICE\Scripting\results\PhotronUSAFCali\spline_max\z_positions.mat';
% If given this folder, the images are extracted from there
ims_folder = 'D:\HoloICE\Holograms\Glass beads\Photron - 2um\20210812\Calibration USAF';

load(path_to_file)

% plot_indices(path_to_file,zs,autof_Ws,'sum',[],legend)   % Filtered

sector_bb = [296, 435, 251, 350];
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

mean_goods = zeros(Np, Ne, Nt,Nim-1);
std_mean_goods = zeros(Np, Ne, Nt,Nim-1);
max_goods = zeros(Np, Ne, Nt,Nim-1);
Ntotal = numel(mean_goods);

msg = 'Processing sweep...';
wbar = waitbar(0, msg);
time_ini = clock;
last_refresh = clock;
refreshing_time = 5;    % seconds

for ki = 1:Npars
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
   clear f1 f2 r Hd win
   
   for kim = 1:Nim
      if exist(ims_folder,'dir')
         image_path = fullfile(ims_folder, fileList(kim).name);
      else
         image_path = fullfile(fileList(kim).folder, fileList(kim).name);
      end
      % Read images
      im = imread(image_path);
      Psi = hologram_reconstruction_at_z(im, z_focus(kim), wavelength, dx);
      aPsi = abs(Psi);
      clear Psi
      aPsi = filter2(h, aPsi);
      
      % SIFT vlFeat
      if kim == 1
         aPsi_min = min(aPsi(:));
         aPsi_max = max(aPsi(:));
      end
      aPsi = single((aPsi-aPsi_min)/(aPsi_max-aPsi_min)*255);
      if kim > 1
         [f,d] = vl_sift(aPsi, 'PeakThresh', peak_thresh, 'edgethresh', edge_thresh);
         [matches, scores] = vl_ubcmatch(d_prev, d, thresh_ubc);
         matchedDistortedXY = [f(1,matches(2,:)); f(2,matches(2,:))]';
         matchedOriginalXY = [f_prev(1,matches(1,:)); f_prev(2,matches(1,:))]';
         [tformBackToOriginal,inlierIdx,err_struct] = my_estimateGeometricTransform2D(matchedDistortedXY,...
            matchedOriginalXY,'similarity',MAX_TOL_ERR,false);
         % show_match_sift(aPsi_prev,aPsi,matches,f_prev,f,inlierIdx,tformBackToOriginal);
         f_prev = f;
         d_prev = d;
         
         mean_goods(kp, ke, kt,kim-1) = err_struct.mean_good;
         std_mean_goods(kp, ke, kt,kim-1) = err_struct.std_good / sqrt(err_struct.n_good);
         max_goods(kp, ke, kt,kim-1) = err_struct.max_good;
      else
         [f_prev,d_prev] = vl_sift(aPsi, 'PeakThresh', peak_thresh, 'edgethresh', edge_thresh);
         % show_im_sift(aPsi, f_prev, [], d_prev)     % Show regions to get features and matching detection.
      end
      aPsi_prev = aPsi;
   end
   
   if etime(clock, last_refresh) > refreshing_time
      last_refresh = clock;
      elapsed_time = etime(clock, time_ini);
      remaining_time = (Npars - ki) * elapsed_time / ki;
      if remaining_time > 120
         msg = sprintf('Processing z sweep... %d minutes remaining.', round(remaining_time/60));
      else
         msg = sprintf('Processing z sweep... %d seconds remaining.', round(remaining_time));
      end
   end
   waitbar(ki/Npars,wbar, msg);
end
close(wbar)
save('results\sweep_pars_SIFT_out.mat','mean_goods','std_mean_goods','max_goods','peaks_th','edges_th','thrs_ubc','path_to_file')
