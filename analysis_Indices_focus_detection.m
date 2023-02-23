% Analysis of z_positions.mat, results from Focus_detection_single_plane.
% Obtains the scaling between different hologram reconstructions.
% The output can be used by focus_to_camera_estimation_USAF
clear variables
close all

dx = 8e-6; % 20e-6;
wavelength = 355e-9; % 532e-9;

path_to_file = 'D:\HoloICE\Holograms\Pollen\12 - Pollen Optronis\Calibration USAF\results\z_positions.mat';
      % 'D:\HoloICE\Holograms\Glass beads\Photron\20210916_Microspheres2um\SeqCalib01\results\z_positions.mat';
ims_path = 'D:\HoloICE\Holograms\Pollen\12 - Pollen Optronis\Calibration USAF';
      % 'D:\HoloICE\Holograms\Glass beads\Photron\20210916_Microspheres2um\SeqCalib01';
PLOT_BOX = false;
do_plot_register = true;
use_SIFT = true;

load(path_to_file)
% z_focus = z_focus(1:end-1);   % Photron 15 um
% fileList = fileList(1:end-1);

% plot_indices(path_to_file,zs,autof_Ws,'sum',[],legend)   % Filtered

% sector_bb = [296, 435, 251, 350];   % Debug with numbers 2 and 3 in first images
% sector_bb = [268, 768, 20, 640];
sector_bb = [150, 700, 1100, 1500];   % Photron 15um
hp_freq = 0.1;

% Create highpass filter
[f1,f2] = freqspace(21,'meshgrid');
r = sqrt(f1.^2 + f2.^2);
Hd = ones(21); 
Hd(r<hp_freq) = 0;
win = fspecial('gaussian',21,2);
win = win ./ max(win(:));
h = fwind2(Hd,win);
clear f1 f2 r Hd win

peak_thresh = 6;  % Parameters for SIFT registration (lower to get more points)
edge_thresh = 7;  % Parameters for SIFT registration (increase to get more points)
thresh_ubc = 2; %  A descriptor D1 is matched to a descriptor D2 only if the distance d(D1,D2) multiplied by THRESH is not greater than the distance of D1 to all other descriptors. The default value of THRESH is 1.5.
MAX_TOL_ERR = 3.5; % Maximum tolerated error when detecting outliers in matched features by comparing with expected transformed location of initial features.

Nim = length(fileList);

moving_to_camera = (z_focus(2) - z_focus(1)) < 0;
sector_bb_xy = [sector_bb(2:-1:1)', sector_bb(4:-1:3)'; [1 1]];

[filepath,~,~] = fileparts(path_to_file);
output_file = fullfile(filepath,'spacings_SIFT.mat');
if exist(output_file,'file')
   load(output_file,'scales_recovered','scales_recovered_SIFT','tforms_forward', 'imBgnd', 'sector_bb')
else
   % Get background to enhance contrast
   BGND_FILTER_FACTOR = 32;
   for kim = 1:Nim
      image_path = fullfile(ims_path, fileList(kim).name);
      im = imread(image_path);
      [imH, imW] = size(im);
      if kim == 1
         imBgnd = zeros(imH, imW, 'single');
      end
      im = single(im);
      imBgnd = imBgnd + imgaussfilt(im, imW/BGND_FILTER_FACTOR);
   end
   imBgnd = imBgnd / Nim;
   
   scales_recovered = zeros(Nim-1,1);
   scales_recovered_SIFT = zeros(Nim-1,1);
   tforms_forward = cell(Nim-1,1);
   % tforms_back = cell(Nim-1,1);
   if moving_to_camera
      sorted_images = 1:Nim;
   else
      sorted_images = Nim:-1:1;
   end
   for kim = 1:Nim   
      kim_ = sorted_images(kim);
      fprintf('Working with image: %d\n', kim_)
      if exist(ims_path,'dir')
         image_path = fullfile(ims_path, fileList(kim_).name);
      else
         image_path = fullfile(fileList(kim_).folder, fileList(kim_).name);
      end
      % Read images
      im = imread(image_path);
      im = single(im) ./ imBgnd;      
      Psi = hologram_reconstruction_at_z(im, z_focus(kim_), wavelength, dx);
      aPsi = abs(Psi);
      clear Psi
      aPsi = filter2(h, aPsi);   

      if PLOT_BOX
         figure, imagesc(aPsi), hold on
         x1 = sector_bb(3); x2 = sector_bb(4);
         y1 = sector_bb(1); y2 = sector_bb(2);
         x_sq = [x1, x2, x2, x1, x1];
         y_sq = [y1, y1, y2, y2, y1];
         plot(x_sq, y_sq, 'b-', 'LineWidth', 3);
         pause
      end

      % SIFT vlFeat
      if kim == 1
         aPsi_min = min(aPsi(:));
         aPsi_max = max(aPsi(:));
      end
      
      aPsi = single((aPsi-aPsi_min)/(aPsi_max-aPsi_min)*255);
      if kim > 1
         % When the sample is moved towards the camera, the magnification
         % get reduced and the image shows a larger area of the sample,
         % with less resolution.
         if use_SIFT
            [f,d] = vl_sift(aPsi, 'PeakThresh', peak_thresh, 'edgethresh', edge_thresh);      
            [matches, scores] = vl_ubcmatch(d_prev, d, thresh_ubc);
            matchedDistortedXY = [f(1,matches(2,:)); f(2,matches(2,:))]';  % The distorted version is the new one (reduced when moving towards camera)
            matchedOriginalXY = [f_prev(1,matches(1,:)); f_prev(2,matches(1,:))]'; % The original version is the previous one
            % Obtain the transformation from the distorted (reduced when moving towards camera) to the previous
            [tformBackToOriginal,inlierIdx,err_struct] = my_estimateGeometricTransform2D(matchedDistortedXY,...
               matchedOriginalXY,'similarity',MAX_TOL_ERR);
%             show_match_sift(aPsi_prev,aPsi,matches,f_prev,f,inlierIdx,tformBackToOriginal);
            % Get scale and rotation
            tformForward = invert(tformBackToOriginal);  % Transformation from the previous to the present
            Tinv = tformForward.T;
            ss = Tinv(2,1);
            sc = Tinv(1,1);
            scales_recovered_SIFT(kim-1) = sqrt(ss*ss + sc*sc);   % <1 when moving towards the camera
            f_prev = f;
            d_prev = d;
         end         
         
         % Scale the window of registration using the full frame tform
         this_sector_bb_xy = tformForward.T' * sector_bb_xy;
         this_sector_bb = round(reshape(flipud(this_sector_bb_xy(1:2, :)),1,[]));
         % Update the initial tform of registration by adding the effect of windowing the fixed image
         this_tform_ini = affine2d([1 0 0; 0 1 0; this_sector_bb([3,1]) 1]);
         tform_ini_forward = this_tform_ini;
         tform_ini_forward.T = tformForward.T / this_tform_ini.T;
         % Window the fixed image considering the scaling
         fixed_patch_Psi = aPsi(this_sector_bb(1):this_sector_bb(2), this_sector_bb(3):this_sector_bb(4));
         
         % Refine by image registration
%          fixed_patch_Psi = aPsi(sector_bb(1):sector_bb(2), sector_bb(3):sector_bb(4)); % Get a piece of the distorted field (reduced when moving towards camera)
%          tform_ini = affine2d([1 0 0; 0 1 0; sector_bb([3,1]) 1]);      
%          % First apply transform to adjust patch extraction and then go back to prev
%          tform_ini_back = tform_ini;
%          if use_SIFT
%             tform_ini_back.T = tform_ini.T * tformBackToOriginal.T;
%          end
%          tform_ini_forward = invert(tform_ini_back);
         imrf = imref2d(size(fixed_patch_Psi));
   %       movingRegistered_ini = imwarp(aPsi_prev,invert(tform_ini),'OutputView',imrf); % Images overlapped but not scaled      
         % Plot initial transform based on SIFT feature matching
   %       movingRegistered_ini_back = imwarp(aPsi_prev,invert(tform_ini_back),'OutputView',imrf); % Images overlapped and scaled
   %       figure, subplot(1,3,1),imshowpair(fixed_patch_Psi,movingRegistered_ini_back)
   %       subplot(1,3,2), imagesc(fixed_patch_Psi), axis image
   %       subplot(1,3,3), imagesc(movingRegistered_ini_back), axis image
         % Image registration
         [optimizer, metric] = imregconfig('multimodal');
   %       optimizer.InitialRadius = 1e-12;    % Default 6.25e-3
         optimizer.Epsilon = 1.5e-4;       % Default 1.5e-6, lower for accuracy
         optimizer.GrowthFactor = 1.01;      % Default 1.05, lower for precision
   %       optimizer.MaximumIterations = 300;  % Default 100
         % Get the transformation for going forward (from the previous to
         % the distorted (reduced when moving towards camera) :
         % inv(backToOriginal) inv(crop)
         tform = imregtform(aPsi_prev,fixed_patch_Psi,'similarity',optimizer,metric,'InitialTransformation',tform_ini_forward);
         if do_plot_register
            movingRegistered = imwarp(aPsi_prev,tform,'OutputView',imrf);
            figure, subplot(1,3,1),imshowpair(fixed_patch_Psi,movingRegistered)
            subplot(1,3,2), imagesc(fixed_patch_Psi), axis image
            subplot(1,3,3), imagesc(movingRegistered), axis image
         end
         ss = tform.T(2,1);
         sc = tform.T(1,1);
         % Only save the part corresponding to the whole image, not the cropped one
         tform.T = tform.T * this_tform_ini.T; 
         tforms_forward{kim-1} = tform;   
         % Save scale
         scales_recovered(kim-1) = sqrt(ss*ss + sc*sc);
      else
         [f_prev,d_prev] = vl_sift(aPsi, 'PeakThresh', peak_thresh, 'edgethresh', edge_thresh); 
%          show_im_sift(aPsi, f_prev, [], d_prev)     % Show regions to get features and matching detection.           
      end
      aPsi_prev = aPsi;
   end
   disp(scales_recovered)
   [filepath,name,ext] = fileparts(path_to_file);
   save(fullfile(filepath,'spacings_SIFT.mat'),'scales_recovered','scales_recovered_SIFT','tforms_forward', 'imBgnd', 'sector_bb')
end

% Getting scales from each image to the following ones (less resolution, moving towards the camera)
if moving_to_camera
   range_reference = 1:Nim-1;
   range_others = 2:Nim;
else
   range_reference = Nim:-1:2;
   range_others = Nim-1:-1:1;
end
scales = zeros(Nim-1); % Rows: 1,...,Nim-1   Cols: 2,...,Nim
tforms = cell(Nim-1);
centers = zeros(Nim-1,Nim-1,2);
[optimizer, metric] = imregconfig('multimodal');
tform_ini = affine2d([1 0 0; 0 1 0; sector_bb([3,1]) 1]);   % Transformation for cropping
optimizer.Epsilon = 1.5e-4;       % Default 1.5e-6, lower for accuracy
optimizer.GrowthFactor = 1.01;      % Default 1.05, lower for precision
metric.NumberOfSpatialSamples = 750;   % Default 500
for kim = 1:Nim-1
   kim_ = range_reference(kim);
   % Get the reference image
   fprintf('Working with image: %d as reference\n', kim_)
   if exist(ims_path,'dir')
      image_path = fullfile(ims_path, fileList(kim_).name);      
   else
      image_path = fullfile(fileList(kim_).folder, fileList(kim_).name);      
   end
   % Read image
   im = imread(image_path);
   im = single(im) ./ imBgnd;
   Psi = hologram_reconstruction_at_z(im, z_focus(kim_), wavelength, dx);
   aPsi = abs(Psi);
   clear Psi
   aPsi_ref = filter2(h, aPsi); % Reference reconstruction (magnified)
   current_forward_T = eye(3); % inv(tform_ini.T);
   
   for klow = kim:Nim-1
      klow_ = range_others(klow);
      fprintf('Scale with respect to image: %d\n', klow_)
      if exist(ims_path,'dir')
         image_path = fullfile(ims_path, fileList(klow_).name);
      else
         image_path = fullfile(fileList(klow_).folder, fileList(klow_).name);
      end
      % Read image
      im = imread(image_path);
      im = single(im) ./ imBgnd;
      Psi = hologram_reconstruction_at_z(im, z_focus(klow_), wavelength, dx);
      aPsi = abs(Psi);
      clear Psi
      aPsi = filter2(h, aPsi);
      
      % Image registration
      tform_ini_forward = tform_ini; % Initialize
      % The saved tform_forward is just the transformation for a single
      % step (from klow-1 to klow), forward means towards the camera.
      % Add previous steps and traslation due to cropping the fixed piece.
      tform_ini_forward.T = tforms_forward{klow}.T * current_forward_T;
      % Scale the window of registration using the full frame tform
      this_sector_bb_xy = tform_ini_forward.T' * sector_bb_xy;      
      this_sector_bb = round(reshape(flipud(this_sector_bb_xy(1:2, :)),1,[]));
      % Update the initial tform of registration by adding the effect of windowing the fixed image
      this_tform_ini = affine2d([1 0 0; 0 1 0; this_sector_bb([3,1]) 1]);
      tform_ini_forward.T = tform_ini_forward.T / this_tform_ini.T;
      % Window the fixed image considering the scaling
      fixed_patch_Psi = aPsi(this_sector_bb(1):this_sector_bb(2), this_sector_bb(3):this_sector_bb(4));            
      % Image registration
      tform = imregtform(aPsi_ref,fixed_patch_Psi,'similarity',optimizer,metric,'InitialTransformation',tform_ini_forward);
      % Save registration tform without cropping
      current_forward_T = tform.T * this_tform_ini.T;
      % Get scales
      ss = tform.T(2,1);
      sc = tform.T(1,1);
      scales(kim,klow) = sqrt(ss*ss + sc*sc);
      % Save tforms
      tforms{kim,klow} = current_forward_T;
      % Get centers. Points that stay in their place. Not sure if they mean
      % something
      centers(kim,klow,:) = (current_forward_T(1:2,1:2)'-eye(2)) \ -current_forward_T(3,1:2)';
      if do_plot_register
         imrf = imref2d(size(fixed_patch_Psi));
         movingRegistered = imwarp(aPsi_ref,tform,'OutputView',imrf);
         figure, subplot(1,3,1),imshowpair(fixed_patch_Psi,movingRegistered)
         subplot(1,3,2), imagesc(fixed_patch_Psi), axis image
         subplot(1,3,3), imagesc(movingRegistered), axis image
      end
   end
end

save(fullfile(filepath,'spacings_all2all.mat'),'scales', 'tforms', 'centers', 'moving_to_camera')

%%% Trying MATLAB matchFeatures
%    % Matching features
%    if kim > 1
%       points = detectSURFFeatures(aPsi);
%       [features,valid_points] = extractFeatures(aPsi,points);
%       indexPairs = matchFeatures(features,features_prev);
%       matchedPoints = valid_points(indexPairs(:,1),:);
%       matchedPoints_prev = valid_points_prev(indexPairs(:,2),:);
%       figure; subplot(1,3,1)
%       showMatchedFeatures(aPsi,aPsi_prev,matchedPoints,matchedPoints_prev);
%       subplot(1,3,2), imagesc(aPsi), hold on, plot(matchedPoints.Location(:,1),matchedPoints.Location(:,2),'o')
%       subplot(1,3,3), imagesc(aPsi_prev), hold on, plot(matchedPoints_prev.Location(:,1),matchedPoints_prev.Location(:,2),'x')
%    else
%       points_prev = detectSURFFeatures(aPsi);
%       [features_prev,valid_points_prev] = extractFeatures(aPsi,points_prev);
%    end
%    aPsi_prev = aPsi;

% angles = zeros(Nim, 1);
% scales = zeros(Nim, 1);
% tform_ini = affine2d([1 0 0; 0 1 0; sector_bb(1:2) 1]);
%    if kim == 1
%       movingPsi = aPsi(sector_bb(1):sector_bb(2), sector_bb(3):sector_bb(4));
%    end
%    % Image registration not working
%    [optimizer, metric] = imregconfig('multimodal');
%    optimizer.InitialRadius = 1e-8;    % Default 6.25e-3
% %    optimizer.Epsilon = 1.5e-4;       % Default 1.5e-6
%    optimizer.GrowthFactor = 1.01;      % Default 1.05
%    optimizer.MaximumIterations = 300;  % Default 100
%    tform = imregtform(movingPsi,aPsi,'similarity',optimizer,metric); % ,'InitialTransformation',tform_ini);
%    movingRegistered = imwarp(movingPsi,tform,'OutputView',imref2d(size(aPsi)));
%    
%    figure
%    imshowpair(aPsi,movingRegistered)
%    
%    u = [0 1]; 
%    v = [0 0]; 
%    [x, y] = transformPointsForward(tform, u, v); 
%    dx = x(2) - x(1); 
%    dy = y(2) - y(1); 
%    angles(kim) = (180/pi) * atan2(dy, dx) ;
%    scales(kim) = 1 / sqrt(dx^2 + dy^2);