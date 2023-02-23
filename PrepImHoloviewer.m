% Script for getting an image ready to process with Holoviewer.
% Crop, decimation, removing bad rows
% It also informs the reconstruction distances and whether they are okay to
% be processed by ASM. Effective wavelength using windows thicknesses

clear variables
close all
mm = 1e-3;
um = 1e-6;
nm = 1e-9;

% image_path = 'C:\Users\j05625pe\Documents\Work\Holograms\Spherical R2L2S3P1\R2L2S3P1_sphe_03.png';
% image_path = 'C:\Users\j05625pe\Documents\Work\Holograms\Spherical R2L2S3P1\20201026\Im03.tif';    
% image_path = 'C:\Users\j05625pe\Documents\Work\Holograms\Spherical R2L2S3P1\20201117\Background00_4pulses.tif';
% image_path = 'C:\Users\j05625pe\Documents\Work\Holograms\Spherical R2L2S3P1\20210215\Background**.tif';     % Use ** to indicate that any with that pattern is good and they will be averaged
image_path = 'C:\Users\j05625pe\Documents\Work\Holograms\Spherical_Levitation\20210329\Drop_White_gauge_notRainX\image_drop.*.tiff';
bgnd_path = 'C:\Users\j05625pe\Documents\Work\Holograms\Spherical_Levitation\20210329\Background\image_drop._down1_aver.png'; % [];

downsampling = 1;
downsampling_allow_auto_correction = false;
crop_factor = 1;
z_autofoc =  [15*mm, 70*mm];  % zeval; % 
L = 78 *mm; % Distance from sensor to focus point.
use_equivalent_wavelength = false;
EQUALIZE = true;
PICK_OLD = true; % true for using old results

dx0 = 8 *um; %  2.2 *um;
dy0 = 8 *um; % 2.2 *um;
wavelength = 355 *nm;

% Fused silica refractive index
% nFS = sqrt(1+.6961663 * wavelength^2/(wavelength^2-.0684043^2) + .4079426 * ...
%    wavelength^2/(wavelength^2-.1162414^2) + .8974794 * wavelength^2/...
%    (wavelength^2-9.896161^2));
nFS = 1.4761; % Fused Silica at 355 nm
nSapphire = 1.793; % Average for extraordinaire and ordinair sapphire at 355 nm
n1 = nSapphire;
t1 = 5 *mm;
n2 = nFS;
t2 = 12 *mm;
% Dealing with windows effects
if use_equivalent_wavelength
   wavelength_equivalent = wavelength * L / (n1*t1+n2*t2 + (L-t1-t2));
else
   t1_eq = t1 / n1;
   t2_eq = t2 / n2;
   L = L - t1 + t1_eq - t2 + t2_eq;
   wavelength_equivalent = wavelength;
   if any(z_autofoc >= L)
      warning("Some z positions from sweep were extracted from the analysis as being greater or equal to the equivalent distance of %.2f mm",L*1e3)
      z_autofoc = z_autofoc(z_autofoc < L);
   end
end

% Reconstruction distances using plane wave
z_recon_autofoc = z_autofoc * L ./ (L-z_autofoc);

% See if there are *
[filepath,name,ext] = fileparts(image_path);
pieces = split(string(name),'*');
images_path = {};
if length(pieces) == 1
   Nfiles = 1;
   images_path = {[name ext]};
else
   % Get the regular expresion of the files
   the_regexp = strjoin(pieces,'.');
   the_regexp = the_regexp + '\' + ext;
   list_files = dir(filepath);
   Nfiles = 0;
   for k=1:length(list_files)
      if strcmp(list_files(k).name, regexp(list_files(k).name,the_regexp,'match'))
         Nfiles = Nfiles + 1;
         images_path{Nfiles} = list_files(k).name;
      end
   end
   root_name = pieces(1);
end

for k = 1:Nfiles
   im_path = fullfile(filepath,images_path{k});
   out_path = [im_path(1:end-4) sprintf('_down%d',downsampling) '.png'];
   if isfile(out_path) && PICK_OLD
      im = double(imread(out_path));
      [Nr, Nc] = size(im); % [Ny, Nx]
   else
      % Read image
      im = imread(im_path);
      % Using only good rows
      [gr_start, gr_end, im_good] = detect_good_row_start_end(im,[],true);
      if ~isempty(im_good)
         im = im_good;
         clear im_good
      end
      im = im(gr_start:gr_end,:);
      [Nr, Nc] = size(im); % [Ny, Nx]
      % Cropping
      Nr_c = floor(Nr * crop_factor); Nc_c = floor(Nc * crop_factor); 
      ir_c = floor(Nr/2- Nr_c/2)+1; ir_c = ir_c:min(ir_c+Nr_c,Nr);
      ic_c = floor(Nc/2- Nc_c/2)+1; ic_c = ic_c:min(ic_c+Nc_c,Nc);
      im = im(ir_c,ic_c);
      [Nr, Nc] = size(im); % [Ny, Nx]
      % Required downsampling to avoid aliasing in all reconstructions
      down_required = max(abs(z_recon_autofoc)) / min(2*Nr*dy0^2,2*Nc*dx0^2) *wavelength_equivalent;   % Zero padding considered
      if down_required > downsampling && downsampling_allow_auto_correction
         downsampling = ceil(down_required);
         warning('Downsampling factor was change to avoid aliasing in reconstruction to value = %.2f.',downsampling);
      end
      % Downsample
      [im, dx, dy, Nc, Nr] = downsample_image(im,downsampling,dx0,dy0);
      % Equalize all images considered by using a background image or the
      % first (or only one).
      if EQUALIZE
         if k == 1
            if isempty(bgnd_path)
               im_bgnd = im;
            else
               im_bgnd = single(imread(bgnd_path));
            end
            n = min(Nr,Nc)/12;
            background = imgaussfilt(im_bgnd,n);
            im = im./background;
            Mim = max(im(:));
            mim = min(im(:));
            eq_gain = 255.99 / (Mim-mim);
            eq_offset = - mim * eq_gain;
            im = im *eq_gain + eq_offset;
         else
            im = im./background *eq_gain + eq_offset;
         end
      end
      imwrite(uint8(im), out_path)
   end
   
   if k==1
      im_aver = im;
   else
      im_aver = im + im_aver;
   end   
end
if Nfiles > 1
   im_aver = im_aver / Nfiles;
   filename_aver = fullfile(filepath,root_name);
   out_path = [char(filename_aver) sprintf('_down%d_aver',downsampling) '.png'];
   imwrite(uint8(im), out_path)
end

fprintf('The reconstruction distances should go from %.2f to %.2f mm.\n',z_recon_autofoc(1)/mm,z_recon_autofoc(2)/mm)
fprintf('The effective wavelength is %.2f nm.\n',wavelength_equivalent/nm)
fprintf('The effective distance is %.2f mm. \n',L/mm)
fprintf('The size of the output image is %d x %d.\n',Nr,Nc)