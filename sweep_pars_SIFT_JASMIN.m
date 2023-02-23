% Cannot make vl_feat work in JASMIN
function sweep_pars_SIFT_JASMIN(index,n_job_array,config_file)
% Analysis of Indices.mat, results from Focus_detection_single_plane.
% z_file = 'C:\Work\Scripting\results\Calibration USAF\first_attempt\z_positions.mat';
% If given this folder, the images are extracted from there
% ims_folder = 'C:\Users\j05625pe\Documents\Work\Holograms\Glass beads\Photron - 2um\20210812\Calibration USAF';

load(config_file,'z_file', 'ims_folder', 'output_file', 'dx', 'wavelength', 'Np', 'Ne', 'Nt', 'Nh', 'bounds_p', 'bounds_e', 'bounds_t', 'bounds_h')
load(z_file, 'z_focus', 'fileList')

% Find ind to work with
index = str2double(index);
n_job_array = str2double(n_job_array);
Npars = Np * Ne * Nt * Nh;
w_ind = 1:Npars;
w_ind = w_ind(ceil(w_ind/Npars*n_job_array) == index);
Ni = length(w_ind);

peaks_th = logspace(log10(bounds_p(1)), log10(bounds_p(2)), Np);
edges_th = logspace(log10(bounds_e(1)), log10(bounds_e(2)), Ne);
thrs_ubc = logspace(log10(bounds_t(1)), log10(bounds_t(2)), Nt);
hp_freqs = logspace(log10(bounds_h(1)), log10(bounds_h(2)), Nh);

MAX_TOL_ERR = 3; % Maximum tolerated error when detecting outliers in matched features by comparing with expected transformed location of initial features.

Nim = length(fileList);

mean_goods = zeros(Ni,Nim-1);
std_mean_goods = zeros(Ni,Nim-1);
max_goods = zeros(Ni,Nim-1);

time_ini = clock;
last_refresh = clock;
refreshing_time = 600;    % seconds

for ki = 1:Ni
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
      image_path = fullfile(ims_folder, fileList(kim).name);
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
         [matches, ~] = vl_ubcmatch(d_prev, d, thresh_ubc);
         matchedDistortedXY = [f(1,matches(2,:)); f(2,matches(2,:))]';
         matchedOriginalXY = [f_prev(1,matches(1,:)); f_prev(2,matches(1,:))]';
         [~,~,err_struct] = my_estimateGeometricTransform2D(matchedDistortedXY,...
            matchedOriginalXY,'similarity',MAX_TOL_ERR,false);
         f_prev = f;
         d_prev = d;
         
         mean_goods(ki,kim-1) = err_struct.mean_good;
         std_mean_goods(ki,kim-1) = err_struct.std_good / sqrt(err_struct.n_good);
         max_goods(ki,kim-1) = err_struct.max_good;
      else
         [f_prev,d_prev] = vl_sift(aPsi, 'PeakThresh', peak_thresh, 'edgethresh', edge_thresh);
      end
      aPsi_prev = aPsi;
   end
   if etime(clock, last_refresh) > refreshing_time
      last_refresh = clock;
      elapsed_time = etime(clock, time_ini);
      remaining_time = (Ni - ki) * elapsed_time / ki;
      if remaining_time > 120
         fprintf('%d minutes remaining.\n', round(remaining_time/60));
      else
         fprintf('%d seconds remaining.\n', round(remaining_time));
      end
   end
end
output_file = sprintf('%s_%05d.mat', output_file, index);
save(output_file,'mean_goods','std_mean_goods','max_goods','w_ind')
fprintf('Done.\n')