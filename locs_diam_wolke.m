clear variables

% This script reads the classification data file or the tree decision file
% to classify particles and extract information about location and diameter
% into a file.

method = 'ASM';
name_folder = '';  % 'whole'; % 'multifocus';   %
name_file = 'local_win_single_ASM.cfg';
type_results = 'SVS';  % 'Photron';   % 'sim';
if strcmp(type_results, 'sim')
   % Configuration of simulation results
   num_particles = 50;
   cfg_file = sprintf('D:\\HoloICE\\MyHolosuite\\Results\\SimsSPIE\\%s\\mixed_size_num_particles_%d\\%s\\%s',method,num_particles,name_folder,name_file);
   nameWolke = sprintf('SimSPIE_%s_%d',method,num_particles);
elseif strcmp(type_results, 'Photron') 
   if strcmp(method,'SAASM')
      method_dir = 'C001H001S0001';
   else
      method_dir = method;
   end
   cfg_file = sprintf('D:\\HoloICE\\MyHolosuite\\Results\\Glass beads\\Photron\\20210916_Microspheres2um\\2um_inv1.5Ms_maxpower_%s\\contrast_seq0018\\%s\\%s',method_dir,name_folder,name_file);
   nameWolke = sprintf('Experimental_%s',method);
else
   cfg_file = 'D:\HoloICE\MyHolosuite\Results\Glass beads\10 - New beads\SVS Vistek 2um\Seq2\contrast_ims_scale1.00';   
   cfg_file = fullfile(cfg_file,name_folder,name_file);
   nameWolke = sprintf('Experimental_%s_%s',type_results,method);
end

[folderWolke,~,~] = fileparts(cfg_file);
no_hist = 1;
save_locs_diam = fullfile(folderWolke, 'locs_diams.mat');
renew_files = true;

% Delete files
if renew_files
   if exist(save_locs_diam,"file")
      delete(save_locs_diam)
   end
   file_paths = {fullfile(folderWolke,"removeRepeatedPart.mat");
      fullfile(folderWolke,nameWolke + "_tsData.mat"); fullfile(folderWolke,nameWolke + "_wolke.mat")};
   for k=1:3
      if exist(file_paths{k},"file")
         delete(file_paths{k})
      end
   end
end

% nube = wolke(cfg_file, nameWolke, folderWolke);
nube = wolke.create(cfg_file, nameWolke, folderWolke, folderWolke, no_hist);     % 3rd input is where the hist files are, 4th input is where the output is saved

% Keep only the particles, not the artifacts if any
if nube.parameter.loadCD
   ids_prtcls = find(nube.pData.catLabel ~= 'Artifact');
else
   ids_prtcls = find(nube.pData.catPredict ~= 'Artifact');
end
fieldnames = fields(nube.pData);
for kf = 1:length(fieldnames)
   nube.pData.(fieldnames{kf}) = nube.pData.(fieldnames{kf})(ids_prtcls);
end
nube.save;

% Bring the dimensions to real positions in case that a spherical
% illumination was used
if nube.cfg.Ldivergent
   % Bring area, majsiz, minsiz, xpos, ypos and zpos back to pixels
   pixelsize   = nube.cfg.dx*nube.cfg.dy;
   nube.pData.area = nube.pData.area / pixelsize;
   nube.pData.majsiz = nube.pData.majsiz / sqrt(pixelsize);
   nube.pData.minsiz = nube.pData.minsiz / sqrt(pixelsize);
   nube.pData.xpos   = nube.pData.xpos / nube.cfg.dx + nube.parameter.Nx/2 + 1;  % in pixels uint
   nube.pData.ypos   = nube.pData.ypos / nube.cfg.dy + nube.parameter.Ny/2 + 1;  % in pixels uint

   % Get the magnification at every position and the real z position
   % starting from the source point
   if nube.cfg.useSAASM
      magnifications = nube.cfg.Ldivergent ./ (nube.cfg.Ldivergent - nube.pData.zpos);
      real_zpos = nube.cfg.Ldivergent - nube.pData.zpos;
   else
      magnifications = (nube.pData.zpos + nube.cfg.Ldivergent) / nube.cfg.Ldivergent;
      real_zpos = nube.cfg.Ldivergent ./ magnifications;
   end

   % Get real positions and features
   pixelsizes = nube.cfg.dx*nube.cfg.dy ./ magnifications.^2;
   nube.pData.area = nube.pData.area .* pixelsizes;
   nube.pData.majsiz = nube.pData.majsiz .* sqrt(pixelsizes);
   nube.pData.minsiz = nube.pData.minsiz .* sqrt(pixelsizes);
   nube.pData.xpos   = (nube.pData.xpos - nube.parameter.Nx/2 -1) ./ magnifications * nube.cfg.dx;
   nube.pData.ypos   = (nube.pData.ypos - nube.parameter.Ny/2 -1) ./ magnifications * nube.cfg.dy;
   % Derived parameters
   nube.pData.eqsiz  = mean([nube.pData.minsiz nube.pData.majsiz],2); % 'Mean of majsiz and minsiz'
   nube.pData.equivdia = (4/pi*nube.pData.area).^(1/2);
   % nube.pData.asprat = nube.pData.majsiz ./ nube.pData.minsiz; % 'Particle Aspect Ratio (majsiz/minsiz)', shouldn't change
   nube.pData.pDiam = nube.pData.pDiam ./ magnifications;
end

locations = [nube.pData.xpos, nube.pData.ypos, real_zpos];
diameters = nube.pData.pDiam;
num_hologram = nube.pData.holonum;
save(save_locs_diam,"locations","diameters","num_hologram");

% Plot histograms
subplot(2,2,1), histogram(real_zpos,20), title('Z position')
subplot(2,2,2), histogram(nube.pData.area,20), title('Area')
subplot(2,2,3), histogram(nube.pData.pDiam,20), title('Diameter (rescaled particle)')
subplot(2,2,4), histogram(nube.pData.asprat,20), title('Aspect ratio (majsiz/minsiz)')

