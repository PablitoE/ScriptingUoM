% Script for getting comparisons of particles
cfg_file_ASM = 'D:\HoloICE\MyHolosuite\Results\Pollen\SVS\Ash_02\contrast_seq0000\ASM\local_win_single.cfg';
cfg_file_SAASM = 'D:\HoloICE\MyHolosuite\Results\Pollen\SVS\Ash_02\contrast_seq0000\SAASM\local_win_single.cfg';
holo_file = 'im_008.png';
save_ims_dir = 'D:\HoloICE\MyHolosuite\Results\Pollen\SVS\Ash_02\contrast_seq0000';
save_ims_name_format = 'compare_%.3f_%d_%d_%s.svg';   % z, x:pixel, y:pixel, method
N_ims = 128;
do_SAASM = true;
Kz_step = 4 * 3;
zs_around = 3;
NumWorkers = 1;

% z_asm [mm], x [um], y[um]
particles_positions_ASM = [% 129.9778, 844, 7928; % 44.2150
                           % 49.3191, 3412, -3546; % 28.4097 Muy chico
                           % 53.5, 6942, 9202; % distancia reducida
                           % 48.6528, -4951, -8942; % 28.1873 OK
                           127.79, 12552, 3468; % distancia modificada
                           % 48.6915, 1544, -4323; % 28.2003 OK
                           % 49.1214, -5664, -587; % 28.3440
                           48.54, 11175, -3883; % 28.15
                           48.54, 11204, -1681; % 28.15
                           % 565.9033, 12077, -6413; % Mucha memoria
                           55.979, -10055, 8456]; % 30.5
p = gcp("nocreate");
if (isempty(p) && NumWorkers > 1) || (~isempty(p) && p.NumWorkers ~= NumWorkers)
   if ~isempty(p)
      delete(p)
   end
   if NumWorkers > 1
      p = parpool('local',NumWorkers);
   end
end
% Read config file and prepare hologram
cfg_ASM = config(cfg_file_ASM);
cfg_ASM.current_holo = holo_file;
if do_SAASM
   cfg_SAASM = config(cfg_file_SAASM);
   cfg_SAASM.current_holo = holo_file;
end
% Get info from particles
n_particles = size(particles_positions_ASM,1);
particles_positions_ASM(:,1) = particles_positions_ASM(:,1) * 1e-3;  % mm to m
particles_positions_ASM(:,2:3) = particles_positions_ASM(:,2:3) * 1e-6; % um to m
dx = cfg_ASM.dx;
dy = cfg_ASM.dy;
xpos = particles_positions_ASM(:,2);
ypos = particles_positions_ASM(:,3);
zs_ASM = particles_positions_ASM(:,1);
zs_SAASM = cfg_ASM.Ldivergent * zs_ASM ./ (cfg_ASM.Ldivergent + zs_ASM);

% Read hologram
preimage = img;
preimage.filter_handle = 'prefilters';
preimage.config_handle = cfg_ASM;
preimage.unregisterListeners;
hologram = preimage.ampFilter(cfg_ASM.full_holo_path);
Nx = size(hologram,2);
Ny = size(hologram,1);

% Prepare z positions
n_ims = zs_around*2+1;
sensor_dim = min(Nx * dx, Ny * dy) / 2;
NA_ASM = sensor_dim ./ sqrt(sensor_dim ^ 2 + zs_ASM .^ 2);
resolution_z_ASM = cfg_ASM.lambda ./ NA_ASM .^ 2 / 2;
all_zs_ASM = zs_ASM + (-zs_around:zs_around) * Kz_step .* resolution_z_ASM;
NA_SAASM = sensor_dim ./ sqrt(sensor_dim ^ 2 + zs_SAASM .^ 2);
resolution_z_SAASM = cfg_ASM.lambda ./ NA_SAASM .^ 2 / 2;
all_zs_SAASM = zs_SAASM + (-zs_around:zs_around) * Kz_step .* resolution_z_SAASM;

% Prepare propagators
prop_ASM = Propagator(); % Propagation object
prop_ASM.should_normalize = true;
prop_ASM.config_handle = cfg_ASM;
prop_ASM.preconstruct(hologram,zs_ASM); % With this hologram   
prop_ASM.FPrepped_root;
if ~cfg_ASM.useBLASM
   prop_ASM.FPrepped_filter;
end
prop_ASM.unregisterListeners;
propcfg_ASM     = prop_ASM.dump_cfg;
propstruct_ASM  = prop_ASM.dump;
propstruct_ASM.config_handle   = propcfg_ASM;
clear prop_ASM

if do_SAASM
   prop_SAASM = Propagator(); % Propagation object
   prop_SAASM.should_normalize = true;
   prop_SAASM.config_handle = cfg_SAASM;
   prop_SAASM.fft_workers = 6;
   prop_SAASM.verbose = false;
   prop_SAASM.preconstruct(hologram,all_zs_SAASM(:)); % With this hologram   
   prop_SAASM.unregisterListeners;
   propcfg_SAASM     = prop_SAASM.dump_cfg;
   propstruct_SAASM  = prop_SAASM.dump;
   propstruct_SAASM.config_handle   = propcfg_SAASM;
   clear prop_SAASM
end
fprintf("Preprocessing completed at %s.\n", datetime(now,'ConvertFrom','datenum'))

% Reconstruct and show
for cnt = 1:n_particles
   % Get range of pixels to show
   xp = round(xpos(cnt) / dx + 1 + Nx/2);
   yp = round(ypos(cnt) / dy + 1 + Ny/2);
   xy_lims = [ceil(xp - N_ims/2), floor(xp + N_ims/2 - 1),ceil(yp - N_ims/2), floor(yp + N_ims/2 - 1)];
   xy_lims([1,3]) = max(1,xy_lims([1,3]));
   xy_lims(2) = min(Nx,xy_lims(2));
   xy_lims(4) = min(Ny,xy_lims(4));

   % Prepare images
   n_ims = zs_around*2+1;
   images_ASM = zeros(N_ims, N_ims, n_ims);
   images_SAASM = zeros(N_ims, N_ims, n_ims);
   this_zs_ASM = all_zs_ASM(cnt, :);
   this_zs_SAASM = all_zs_SAASM(cnt, :);
   if ~isempty(gcp("nocreate"))
      parfor cnt_focus = 1:n_ims
         % Do the ASM reconstruction/propagation
         prop_ASM = Propagator(propstruct_ASM);
         [slice_ASM, ~]   = prop_ASM.slice(this_zs_ASM(cnt_focus));
         % Do the ASM reconstruction/propagation
         if do_SAASM
            prop_SAASM = Propagator(propstruct_SAASM);
            [slice_SAASM, ~]   = prop_SAASM.slice(this_zs_SAASM(cnt_focus));
         end

         % Get small images
         images_ASM(:,:,cnt_focus) = abs(slice_ASM(xy_lims(3):xy_lims(4),xy_lims(1):xy_lims(2)));
         if do_SAASM
            images_SAASM(:,:,cnt_focus) = abs(slice_SAASM(xy_lims(3):xy_lims(4),xy_lims(1):xy_lims(2)));
         end
      end
   else
      for cnt_focus = 1:n_ims
         % Do the ASM reconstruction/propagation
         prop_ASM = Propagator(propstruct_ASM);
         [slice_ASM, ~]   = prop_ASM.slice(this_zs_ASM(cnt_focus));
         % Do the ASM reconstruction/propagation
         if do_SAASM
            prop_SAASM = Propagator(propstruct_SAASM);
            [slice_SAASM, ~]   = prop_SAASM.slice(this_zs_SAASM(cnt_focus));
         end

         % Get small images
         images_ASM(:,:,cnt_focus) = abs(slice_ASM(xy_lims(3):xy_lims(4),xy_lims(1):xy_lims(2)));
         if do_SAASM
            images_SAASM(:,:,cnt_focus) = abs(slice_SAASM(xy_lims(3):xy_lims(4),xy_lims(1):xy_lims(2)));
         end
      end
   end

   for cnt_focus = 1:n_ims
      subplot(2,zs_around+1,cnt_focus)
      imagesc(images_ASM(:,:,cnt_focus)), axis image, colormap(bone(256))
      title(sprintf("ASM z position = %.3f mm", this_zs_ASM(cnt_focus)*1e3))
   end
   if cnt == 1
      set(gcf,'OuterPosition',[56, 123, 1725, 876])
   end
   im_name = sprintf(save_ims_name_format,zs_ASM(cnt)*1e3,xp,yp,'ASM');
   saveas(gcf,fullfile(save_ims_dir,im_name))

   if do_SAASM
      for cnt_focus = 1:n_ims
         subplot(2,zs_around+1,cnt_focus)
         imagesc(images_SAASM(:,:,cnt_focus)), axis image, colormap(bone(256))
         title(sprintf("SAASM z position = %.3f mm", this_zs_SAASM(cnt_focus)*1e3))
      end
      im_name = sprintf(save_ims_name_format,zs_SAASM(cnt)*1e3,xp,yp,'SAASM');
      saveas(gcf,fullfile(save_ims_dir,im_name))
   end
end