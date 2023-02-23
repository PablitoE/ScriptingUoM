% Estimate distance from camera to focus point by taking 2D Fourier peak of
% the samples with grid of points at different zs (known difference)
% The magnification can be found and therefore the other distances.
% A distances text file must be provided with the positions of the grid
% with respect to any position in z on the side of the source. Distances
% should be in mm and 1 per line.

clear variables
close all
mm = 1e-3;
um = 1e-6;
nm = 1e-9;

do_plot = false;

% Camera info
dpix = 8*um; % 2.2*um;    % 20*um;  % 
downsampling = 1; % Allows faster processing
minimum_reps_pattern = 8;  % The images must show at least some periods of the pattern given by the gridded sample

% Grid info
real_spacing = 125*um;

% Difference of distance of target-sensor between images
distances_textfile = 'D:\HoloICE\Holograms\Cold room tests\Nanodispenser Needle\2022_05_16 calibration.txt';
table = readtable(distances_textfile);
distances = table2array(table) * mm;
exclude_last = 11;

% Images of a gridded sample with different length sensor-to-sample
path_ims = 'D:\HoloICE\Holograms\Cold room tests\Nanodispenser Needle\2022_05_16 calibration';

% Substrate information
% Material : soda lime glass
% windows = {'Source: Sapphire_5mm'; 'Sensor: Sapphire_1mm'};
windows = {'Source: FS_5mm'; 'Sensor: FS_12mm'};
nSLG = 1.544;
thSLG = 1.5 *mm;

% PRepare reading images
list_dir = dir(path_ims);
file_format = 'tiff';
valid_file = zeros(length(list_dir),1,'logical');
if isempty(list_dir)
   error('No png files. Try : mogrify -format png *.*')
end
for k=1:length(list_dir)
   valid_file(k) = length(list_dir(k).name) > 4 && strcmp(list_dir(k).name(end-length(file_format)+1:end),file_format);
end
list_dir = list_dir(valid_file);
n_ims = length(list_dir);
names_ims = cell(n_ims,1);
freqs = zeros(n_ims,1);
spacings = NaN*ones(1,n_ims);    % Spatial spacing of features
err = zeros(n_ims,1);
% Try to recover previous results of frequency estimation
info_file = fullfile(path_ims,'windows_distances_spacings.mat');
if exist(info_file, 'file')
   load(info_file)
   [distances, ind_sort] = sort(distances);
   spacings = spacings(ind_sort);
end
work_with = find(isnan(spacings));
for k = work_with
   file = list_dir(k);
   fullpath_file = fullfile(path_ims,file.name);
   names_ims{k} = file.name;

   % Get frequency
   im = imread(fullpath_file);
   [im, dx, dy, Nc, Nr] = downsample_image(im,downsampling,dpix,dpix);
%          % Trying autocorrelation with a subimage... too smooth to find a
%          % good solution
%          Nrsub = round(Nr/12);
%          Ncsub = round(Nc/12);
%          subim = im(1:Nrsub, 1:Ncsub);
%          autoc = xcorr2(im, subim);
   guard = 2 / (Nc / minimum_reps_pattern);
   % (image, minimum frequency, number of maximae to consider (taking the closer to origin), amount of pixels to approx with parabola)
   [freqs(k), err(k)] = get_main_frequency(im,guard,20,2,do_plot,false);
   fprintf('Estimation of period in file: %d\n', k)
   spacings(k) = 2*dx/freqs(k);
end

save(info_file, 'windows', 'spacings', 'distances', "dpix")

if exclude_last > 0
   spacings = spacings(1:end-exclude_last);
   spacings = spacings(:);
   distances = distances(1:end-exclude_last);
end

% Distance of reconstruction from 2019 - Resolution and sampling analysis in digital in-line holography with spherical wave illumination
% L: source–sensor distance
% z : object–sensor distance
% d : object-source distance
% m : L/(L-z) magnification : L/d

% Model with unknown L, d0, real_spacing
% real_spacing * m = real_spacing * L / (d0 + dn) = spacing
% min sum(f(x)^2), x = [rs, L, d0]
model1 = @(x) x(1)*x(2)./(x(3) + distances) - spacings;
x0 = [real_spacing, 50*mm, 10*mm];
lb = [100*um, 20*mm, 0];
ub = [150*um, 250*mm, 220*mm];
options = optimoptions(@lsqnonlin, 'Display','iter', 'OptimalityTolerance', 1e-15, 'FunctionTolerance', 1e-13, 'StepTolerance', 1e-8);
x1 = lsqnonlin(model1,x0,lb,ub, options);

% Model with unknown L, d0, but known real spacing
% real_spacing * m = real_spacing * L / (d0 + dn) = spacing
% min sum(f(x)^2), x = [L, d0]
model2 = @(x) real_spacing*x(1)./(x(2) + distances) - spacings;
x0 = [50*mm, 10*mm];
lb = [20*mm, -100*mm];
ub = [250*mm, 220*mm];
x2 = lsqnonlin(model2,x0,lb,ub,options);

% Using l1-norm to exclude outliers
l1_model = @(x) sum(abs(model2(x)));
[x_l1,fval_l1,exitflag_l1,output_l1] = fsolve(l1_model, x0);

L = x2(1); % Model 1 is not good for estimating L

plot(model1(x1)), hold on
plot(model2(x2))
plot(model2(x_l1))
title('Errors')

figure, plot(x1(3) + distances, x1(2)./(x1(3) + distances)), hold on
plot(x2(2) + distances,x2(1)./(x2(2) + distances))
plot(x_l1(2) + distances,x_l1(1)./(x_l1(2) + distances))
title("Magnification")

% Model vs measurements
figure, plot(distances, x1(2)./(x1(3) + distances) * x1(1)), hold on
plot(distances, x2(1)./(x2(2) + distances) * real_spacing)
plot(distances, x_l1(1)./(x_l1(2) + distances) * real_spacing)
plot(distances, spacings, 'o')
title("Model vs measurements (spacing of dots), spacing : 125 um")

n_SLG = 1.544;
t_grid = 1.5*mm;
L_grid_out = L + t_grid - t_grid / n_SLG;

% These are equivalent distances. Input values to Single_ASM or PrepImHoloviewer should be
nFS = 1.4761; % Fused Silica at 355 nm
nSap = 1.793; % Sapphire
if strcmp(windows{1},'Source: FS_5mm')
   n1 = nFS;
   t1 = 5 *mm;
elseif strcmp(windows{1},'Source: Sapphire_5mm')
   n1 = nSap;
   t1 = 5 *mm;
elseif strcmp(windows{1},'Source: Sapphire_2.3mm')
   n1 = nSap;
   t1 = 2.3 *mm;
else 
   error('Unknown source side window.');
end
if strcmp(windows{2},'Sensor: FS_12mm')
   n2 = nFS;
   t2 = 12 *mm;
elseif strcmp(windows{2},'Sensor: Sapphire_1mm')
   n2 = nSap;
   t2 = 1 *mm;
elseif strcmp(windows{2},'Sensor: Sapphire_2.3mm')
   n2 = nSap;
   t2 = 2.3 *mm;
else 
   error('Unknown sensor side window.');
end
t1_eq = t1 / n1;
t2_eq = t2 / n2;
L_real = L +t1-t1_eq+t2-t2_eq;
% t3_eq = thSLG / nSLG;
% L_real = L +t1-t1_eq+t2-t2_eq+thSLG-t3_eq;

fprintf('The effective distance from the source to the sensor is %.2f mm (%.2f mm with the grid in).\n',L_grid_out/mm,L/mm)
fprintf('The actual distance from the source to the sensor is %.2f mm.\n',L_real/mm)
