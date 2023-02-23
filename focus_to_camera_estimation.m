% Estimate distance from camera to focus point by taking 2D Fourier peak of
% the samples with grid of points at different zs (known difference)
% The magnification can be found and therefore the other distances.

clear variables
close all
mm = 1e-3;
um = 1e-6;
nm = 1e-9;

% Images of a gridded sample with different length sensor-to-sample
% The distance sample-to-focus did not change.
path_ims = '
path_im1 = 'C:\Users\j05625pe\Documents\Work\Holograms\Spherical R2L2S3P1\20201026\Im03.tif'; 
path_im2 = 'C:\Users\j05625pe\Documents\Work\Holograms\Spherical R2L2S3P1\20201026\Im02.tif'; % +1/2"
gr_start1 = 79; % Good rows start : 20201026\Im03:79 20201026\Im02:75
gr_start2 = 75;

% Difference of total length
dL = 12.7 *mm;

dpix = 2.2*um;
downsampling = 4;
real_spacing = 125*um;

im = imread(path_im1);
im = im(gr_start1:end,:);
[im, dx, dy, Nc, Nr] = downsample_image(im,downsampling,dpix,dpix);
freq1 = get_main_frequency(im);
spacing1 = 2*dx/freq1;

im = imread(path_im2);
im = im(gr_start2:end,:);
[im, dx, dy, Nc, Nr] = downsample_image(im,downsampling,dpix,dpix);
freq2 = get_main_frequency(im);
spacing2 = 2*dx/freq2;

fprintf('The magnification factor of changing the position is %.2f.\n', freq1/freq2)

% Distance of reconstruction from 2019 - Resolution and sampling analysis in digital in-line holography with spherical wave illumination
% L: source–sensor distance
% z : object–sensor distance
% d : object-source distance
% m : L/(L-z) magnification : L/d

m1 = spacing1/real_spacing;
m2 = spacing2/real_spacing;

% System of equations
% m1 * d - L = 0
% m2 * d - L = dL

d = -dL/(m1-m2);
L = -dL*m1/(m1-m2);
z = L-d;

% These are equivalent distances. Input values to Single_ASM or PrepImHoloviewer should be
nFS = 1.4761; % Fused Silica at 355 nm
n1 = nFS;
t1 = 2.3 *mm;
n2 = nFS;
t2 = 11.8 *mm;
t1_eq = t1 / n1;
t2_eq = t2 / n2;
L_real = L +t1-t1_eq+t2-t2_eq;
z_real = z +t2-t2_eq;

fprintf('The actual distance from the source to the sensor is %.2f mm.\n',L_real/mm)
fprintf('The actual distance from the sensor to the object is %.2f mm, but the reconstruction distance at the original wavelength is %.2f mm.\n',z_real/mm, z/mm)