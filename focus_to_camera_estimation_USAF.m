% Estimate distance from camera to focus point by taking holograms of a
% USAF target at different positions. % The distance sample-to-focus did
% not change between holograms. The focus plane is found by
% reconstructing the holograms at several positions. The scaling between
% reconstructions are found by SIFT feature matching.

clear variables
close all
mm = 1e-3;
um = 1e-6;
nm = 1e-9;

do_plot = false;

% Camera info
dpix = 8*um;   % 20*um;
Nx = 1696; % 1024;
Ny = 1708; % 1024;
downsampling = 1; % Allows faster processing

% File with relevant estimated values that were outputs from Focus_detections_single_plane.m and analysis_Indices_focus_detection.m 
path_to_files = 'D:\HoloICE\Holograms\Pollen\12 - Pollen Optronis\Calibration USAF\results';
      % 'D:\HoloICE\Holograms\Glass beads\Photron\20210916_Microspheres2um\SeqCalib00\results';
reconstructed_distances_file = 'z_positions.mat';
estimated_spacings_file = 'spacings_all2all.mat'; % 'spacings_SIFT.mat';

% Load variables: z_focus [m], scales_recovered
load(fullfile(path_to_files, reconstructed_distances_file))
load(fullfile(path_to_files, estimated_spacings_file))
% z_focus = z_focus(1:end-1);

if exist('scales_recovered','var')
   z_ini = z_focus(1:end-1);
   z_end = z_focus(2:end);
elseif exist('scales','var')
   ind_ut = triu(ones(size(scales), 'logical'));
   scales_recovered = scales(ind_ut(:));
   z_ini = z_focus(1:end-1);
   z_end = z_focus(2:end);
   [z_end,z_ini] = meshgrid(z_end, z_ini);
   z_ini = z_ini(ind_ut(:));
   z_end = z_end(ind_ut(:));
   cx = cellfun(@(x)get_cx(x), tforms);
   cx = cx(ind_ut(:)) * dpix;
   cy = cellfun(@(x)get_cy(x), tforms);
   cy = cy(ind_ut(:)) * dpix;
else
   error('Unexpected variables loaded from %s.', estimated_spacings_file)
end

% Distance of reconstruction from 2019 - Resolution and sampling analysis in digital in-line holography with spherical wave illumination
% L: source–sensor distance
% z : object–sensor distance
% d : object-source distance
% m : L/(L-z) magnification : L/d
% zr : reconstruction distance for z : z = zr L / (L + zr) or zr = z L / (L - z)

% Model with unknown L
% mag_2 / mag_1 = scale12
% (L-z1) / (L-z2) = scale12
% (1-scale12) L = scale12 zr1 - zr2
% min sum(f(L)^2) : LSE
a = 1 - scales_recovered;
b = scales_recovered.*z_ini - z_end;
L = a'*b/(a'*a);

errors = b./a - L;
RMSE = sqrt(mean(errors.^2));
plot(errors)
title(sprintf('Errors (RMSE = %f)', RMSE))

% m = L / (L - zr/m) = Lm / (Lm - zr)
% Lm - zr = L => m = (L + zr) / L
magnifications = (L+z_focus)/L;
figure, plot(magnifications)
title("Magnification")

% Model vs measurements
% scale12 = (L + zr2) / (L + zr1)
modelled_scalings = (L + z_end) ./ (L + z_ini);
figure, plot(scales_recovered,'x'), hold on
plot(modelled_scalings, 'o')
title("Model vs measurements (scalings)")
legend('Measurements of scales', 'Scales by model')

% Error(L)
Ls = linspace(0.03, 0.2, 200);
errors_L = b./a - Ls;
RMSEs = sqrt(mean(errors_L.^2));
RMSEs_axb = sqrt(mean((a*Ls - b).^2));
figure, plot(Ls, RMSEs, Ls, RMSEs_axb)
title('Errors: RMSE(L)')

% If the reference is not in the center of the hologram then the estimation
% of the point is not directly [x,y,z]/M with [x,y] centered.
% Model (from 9.5.1 Goodman, Introduction to Fourier Optics):
% - plane of sample o: [xo,yo,zo]
% - plane of reference source r: [xr,yr,zr] (zr = L in previous model)
% - plane of reconstruction i: [xi,yi,zi], using reconstruction wave zp=inf
% - zi = M zo; xi = M (xo + xr) (same for y); M = (zr + zi)/zr
% From registration of reconstructions of same sample at zi and zi' (known):
% - xi' = s xi + cx (same for y), for all xo
% - M' (xo + xr) = s [M (xo + xr)] + cx, for all xo, particularly xo = 0:
% - M' xr = sM xr + cx => [(1-s) zr + zi' - s zi] xr - cx zr = 0
% min sum(f(v)^2), v = [xr, yr, zr]
s = scales_recovered; 
model = @(v) [((1-s)*v(3,:) + z_end - s .*z_ini) .* v(1,:) - cx * v(3,:); ((1-s)*v(3,:) + z_end - s .*z_ini) .* v(2,:) - cy * v(3,:)];
x0 = [(Nx+1)/2*dpix, (Ny+1)/2*dpix, -300*mm]';
lb = [2/5*Nx*dpix, 2/5*Ny*dpix, -600*mm]; % [dpix, dpix, 20*mm];
ub = [3/5*Nx*dpix, 3/5*Ny*dpix, 150*mm];
options = optimoptions(@lsqnonlin, 'Display','iter', 'OptimalityTolerance', 1e-12, 'FunctionTolerance', 1e-11);
% options = optimoptions(@lsqnonlin, 'Algorithm', 'levenberg-marquardt', 'Display','iter', 'OptimalityTolerance', 1e-12, 'ScaleProblem', 'jacobian', 'FunctionTolerance', 1e-11);
v_out = lsqnonlin(model,x0,lb,ub,options);
fprintf('Model using offset of reference point source:\n[x,y,z] = [%.2f um, %.2f um, %.2f mm]\n', v_out(1)/um, v_out(2)/um, v_out(3)/mm)

% Take a look at the errors
Nsim_x = 1000;
Nsim_z = 1000;
% dim_xy = linspace(2/5*Nx*dpix,3/5*Nx*dpix,Nsim_x);
% dim_z = linspace(-500*mm, 0*mm, Nsim_z);
dim_xy = linspace(-Nx*dpix,Nx*dpix,Nsim_x);
dim_z = linspace(-500*mm, 500*mm, Nsim_z);
[all_xy, all_z] = meshgrid(dim_xy,dim_z);
all_xy = all_xy(:)';
all_z = all_z(:)';
error = model([all_xy; (Ny+1)/2*dpix*ones(1,Nsim_x*Nsim_z); all_z]);
error = sqrt(mean(error.^2));
error = reshape(error, Nsim_x, Nsim_z);
figure, imagesc(error)

% Analysis considering refractive index of materials
% Bulk material information of USAF target : soda lime glass
nSLG = 1.5261; % @532nm, 1.544; @355nm
thSLG = 1.5 *mm;
t_eq_USAF = thSLG / nSLG;

% Removing USAF target from setup effect on L
L_without_USAF = L + thSLG - t_eq_USAF;
fprintf('The effective distance from the source to the sensor (with the windows, but without the USAF target) is %.2f mm.\n',L_without_USAF/mm);

% These are equivalent distances. Input values to Single_ASM or PrepImHoloviewer should be
% Windows Edmund Optics #39233 Window Saph 25.4 Dia MgF2 Ctd, 2.3 mm thickness
nFS = 1.4761; % Fused Silica at 355 nm
nSap = 1.7717; % Sapphire@532nm, 1.793; % Sapphire@355nm
n1 = nFS; % nSap;
t1 = 5*mm; % 2.3 *mm;
n2 = nFS; % nSap; % 
t2 = 12*mm; % FS % 2.3 *mm; % 
t1_eq = t1 / n1;
t2_eq = t2 / n2;
L_real = L_without_USAF +t1-t1_eq+t2-t2_eq;

fprintf('The actual distance from the source to the sensor is %.2f mm.\n',L_real/mm)

pattern_side_sensor = true; % TODO: think about how this affects the real distance
measured_positions = [22.46 22.76 23.38 23.89 24.39 25.00 25.58 26.00];
z_positions = z_focus .* magnifications;
delta_positions = abs(diff(z_positions)')
delta_measured = diff(measured_positions)

function cx = get_cx(tform)
   if isempty(tform)
      cx = 0; 
   else
      cx = tform(3,1);
   end
end

function cy = get_cy(tform)
   if isempty(tform)
      cy = 0; 
   else
      cy = tform(3,2);
   end
end