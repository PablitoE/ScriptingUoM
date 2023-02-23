% Estimate distance from camera to focus point by taking holograms of a
% USAF target at different positions. % The distance sample-to-focus did
% not change between holograms. The focus plane is found by
% reconstructing the holograms at several positions. The position of the
% sample is measured with respect to an unknown reference.
%
% TODO: verify with new measurements. 

clear variables
close all
mm = 1e-3;
um = 1e-6;
nm = 1e-9;

do_plot = false;

% Camera info
dpix = 20*um;
downsampling = 1; % Allows faster processing

% File with relevant estimated values that were outputs from Focus_detections_single_plane.m and analysis_Indices_focus_detection.m 
path_to_files = 'C:\Work\Scripting\results\Calibration USAF\weigthed_average';
reconstructed_distances_file = 'z_positions.mat';
estimated_spacings_file = 'spacings_SIFT.mat';

% Load variables: z_focus [m], scales_recovered
load(fullfile(path_to_files, reconstructed_distances_file))
load(fullfile(path_to_files, estimated_spacings_file))

% Measured positions of USAF target with respect to a unknown reference
measured_positions = [22.46 22.76 23.38 23.89 24.39 25.00 25.58 26.00];
measured_positions = (measured_positions - measured_positions(1))*mm;
if diff(measured_positions(1:2)) < 0 && diff(z_focus(1:2)) > 0
   measured_positions = - measured_positions;
end

% Distance of reconstruction from 2019 - Resolution and sampling analysis in digital in-line holography with spherical wave illumination
% L: source–sensor distance
% z : object–sensor distance
% d : object-source distance
% m : L/(L-z) magnification : L/d
% zr : reconstruction distance for z : z = zr L / (L + zr) or zr = z L / (L - z)
% dz1 : distance with respect to first z : dz1 = zrL/(L+zr) - zr1L/(L+zr1)

% Model with unknown L
% zrL/(L+zr) - zr1L/(L+zr1) = zmi
% min sum(f(x)^2), x = [rs, L, d0]
zr = z_focus(2:end); zr = zr(:);
zr1 = z_focus(1);
zm = measured_positions(2:end); zm = zm(:);
model = @(x) zr*x./(x+zr) - zr1*x/(x+zr1) - zm;
x0 = 60*mm;
lb = 20*mm;
ub = 250*mm;
options = optimoptions(@lsqnonlin, 'Display','iter', 'OptimalityTolerance', 1e-9, 'FunctionTolerance', 1e-11);
L = lsqnonlin(model,x0,lb,ub,options);

syms zr_ L_ zr1_ zm_
eq = zr_*L_/(L_+zr_) - zr1_*L_/(zr1_+L_) - zm_ == 0;
Lfun = symfun(solve(eq,L_),[zr_,zr1_,zm_]);
Ls = Lfun(zr,zr1,zm); % First cell corresponds to first solution, different inputs

% The error between the result from a single (zri,zr1,zmi) input and the
% lsqnonlin result. Everything from here on must be modified as it was
% copied from focus_to_camera_estimation_USAF

errors = b./a - L;
MSE = sqrt(mean(errors.^2));
plot(errors)
title(sprintf('Errors (MSE = %f)', MSE))

% m = L / (L - zr/m) = Lm / (Lm - zr)
% Lm - zr = L => m = (L + zr) / L
magnifications = (L+z_focus)/L;
figure, plot(magnifications)
title("Magnification")

% Model vs measurements
% scale12 = (L + zr2) / (L + zr1)
modelled_scalings = (L + z_focus(2:end)) ./ (L + z_focus(1:end-1));
figure, plot(scales_recovered,'x'), hold on
plot(modelled_scalings, 'o')
title("Model vs measurements (scalings)")
legend('Measurements of scales', 'Scales by model')

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
n1 = nSap;
t1 = 2.3 *mm;
n2 = nSap; % nFS;
t2 = 2.3 *mm; % 12 FS
t1_eq = t1 / n1;
t2_eq = t2 / n2;
L_real = L_without_USAF +t1-t1_eq+t2-t2_eq;

fprintf('The actual distance from the source to the sensor is %.2f mm.\n',L_real/mm)

pattern_side_sensor = true; % TODO: think about how this affects the real distance
z_positions = z_focus .* magnifications;
delta_positions = abs(diff(z_positions)')
delta_measured = diff(measured_positions)