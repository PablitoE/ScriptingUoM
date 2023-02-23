% Function for getting a set of images ready to process with HoloSuite.
% Crop, decimation, removing bad rows
% It also informs the reconstruction distances and whether they are okay to
% be processed by ASM.
% A background image could be used as reference to get clean holograms.

function PrepSeqHoloSuite(image_path,bgnd_path)

mm = 1e-3;
nm = 1e-9;

% Use ** to indicate that any with that pattern is good and they will be averaged
% image_path = 'C:\Users\j05625pe\Documents\Work\Holograms\Spherical_Levitation\20210429\Seq1\image_drop.*.tiff';
% bgnd_path = 'C:\Users\j05625pe\Documents\Work\Holograms\Spherical_Levitation\20210429\Background\image_drop.*.tiff'; % [];

downsampling = 1;
crop_factor = 1;
z_autofoc =  [30*mm, 80*mm];  % zeval; % 
L = 93 *mm; % Distance from sensor to focus point.
use_equivalent_wavelength = false;
EQUALIZE = true; % Self equalize when analysing images without available background images, otherwise, use background image
PICK_OLD = true; % true for using old results

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

if exist('bgnd_path','var') && ~isempty(bgnd_path)
   % Work with background
   [bgnd_filepath, bgnd_images_names, bgnd_root_name] = get_files_to_process(bgnd_path);
   
   % Get an average image (if does not exist) and save (or load)
   bgnd_aver = get_average_image(bgnd_filepath,bgnd_images_names,bgnd_root_name,downsampling,crop_factor,PICK_OLD);
else
   bgnd_aver = [];
end

% Work with sequence
[images_path, images_names, ~] = get_files_to_process(image_path);

[Nr,Nc] = get_sequence_to_png(images_path,images_names,bgnd_aver,downsampling,crop_factor,PICK_OLD,EQUALIZE);

fprintf('The reconstruction distances should go from %.2f to %.2f mm.\n',z_recon_autofoc(1)/mm,z_recon_autofoc(2)/mm)
fprintf('The effective wavelength is %.2f nm.\n',wavelength_equivalent/nm)
fprintf('The effective distance is %.2f mm. \n',L/mm)
fprintf('The size of the output images is %d x %d.\n',Nr,Nc)

% Move raw images to Raw directory
raw_folder = fullfile(images_path, 'Raw');
if ~isfolder(raw_folder)
   mkdir(raw_folder);
end
for k = 1:numel(images_names)
   movefile(fullfile(images_path,images_names{k}),fullfile(raw_folder,images_names{k}))
end
end

function [filepath, images_names, root_name] = get_files_to_process(image_path)
% See if there are *
[filepath,name,ext] = fileparts(image_path);
pieces = split(string(name),'*');
images_names = {};
if length(pieces) == 1
   images_names = {[name ext]};
else
   % Get the regular expresion of the files
   the_regexp = strjoin(pieces,'.');
   the_regexp = the_regexp + '\' + ext;
   list_files = dir(filepath);
   Nfiles = 0;
   for k=1:length(list_files)
      if strcmp(list_files(k).name, regexp(list_files(k).name,the_regexp,'match'))
         Nfiles = Nfiles + 1;
         images_names{Nfiles} = list_files(k).name;
      end
   end
   root_name = pieces(1);
end
end

function average_image = get_average_image(filepath,images_names,root_name,downsampling,crop_factor,PICK_OLD)
filename_aver = fullfile(filepath,root_name);
aver_out_path = [char(filename_aver) sprintf('_down%d_aver',downsampling) '.png'];
if isfile(aver_out_path) && PICK_OLD
   average_image = double(imread(aver_out_path));
   return
end
Nfiles = length(images_names);
for k = 1:Nfiles
   im_path = fullfile(filepath,images_names{k});
   out_path = [im_path(1:end-4) sprintf('_down%d',downsampling) '.png'];
   if isfile(out_path) && PICK_OLD
      im = double(imread(out_path));
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
      % Downsample
      [im, ~, ~, ~, ~] = downsample_image(im,downsampling);
      % Save
      imwrite(uint8(im), out_path)
   end   
   if k==1
      average_image = im;
   else
      average_image = im + average_image;
   end   
end
if Nfiles > 1
   average_image = average_image / Nfiles;
%    imwrite(uint8(255.99*mat2gray(average_image)), aver_out_path)
   imwrite(uint8(255.99*average_image/max(average_image(:))), aver_out_path)
end
end

function [Nr,Nc] = get_sequence_to_png(filepath,images_names,im_bgnd,downsampling,crop_factor,PICK_OLD,EQUALIZE)
Nfiles = length(images_names);
if ~isempty(im_bgnd)
   [Nr_bgnd, Nc_bgnd] = size(im_bgnd);
end
for k = 1:Nfiles
   im_path = fullfile(filepath,images_names{k});
   out_path = [im_path(1:end-4) sprintf('_down%d',downsampling) '.png'];
   if isfile(out_path) && PICK_OLD
      [Nr, Nc] = size(im_bgnd);
      continue
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
      % Downsample
      [im, ~, ~, Nc, Nr] = downsample_image(im,downsampling);
      % Equalize all images considered by using a background image or the
      % first (or only one).
      if EQUALIZE
         if k == 1
            if isempty(im_bgnd)
               im_bgnd = im;
               n = min(Nr,Nc)/12;
               im_bgnd = imgaussfilt(im_bgnd,n);
            elseif Nr_bgnd ~= Nr || Nc_bgnd ~= Nc
               error('The size of the images is different from the size of background.')
            end  
            im_bgnd(im_bgnd==0) = 1;
            im = im./im_bgnd;
%             Mim = max(im(:));
%             mim = min(im(:));
            qs = quantile(im(:),[0.001, 0.999]);
            mim = qs(1); Mim = qs(2);
            mim = 0; % HARDCODED : 0 amplitude acceptable
            eq_gain = 255.99 / (Mim-mim);
            eq_offset = - mim * eq_gain;
            im = im *eq_gain + eq_offset;
         else
            im = im./im_bgnd *eq_gain + eq_offset;
         end
      end
      imwrite(uint8(im), out_path)
   end  
end
end