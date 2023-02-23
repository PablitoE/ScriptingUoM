clear variables

cfg_file = 'D:\HoloICE\MyHolosuite\Results\SimsSPIE\ASM\mixed_size_num_particles_20\whole\local_win.cfg';
% 'D:\HoloICE\MyHolosuite\Results\Glass beads\10 - New beads\Optronis 2um\Seq3\mode_0\contrast_seq0001\whole\local_win.cfg';
nameWolke = 'SimSPIE_ASM_20';
folderWolke = 'D:\HoloICE\MyHolosuite\Results\SimsSPIE\ASM\mixed_size_num_particles_20\whole';
% 'D:\HoloICE\MyHolosuite\Results\Glass beads\10 - New beads\Optronis 2um\Seq3\mode_0\contrast_seq0001\whole';
no_hist = 1;
save_nube = fullfile(folderWolke, 'nube.mat');

if exist("save_nube", "file")
   load(save_nube)
else
   % nube = wolke(cfg_file, nameWolke, folderWolke);
   nube = wolke.create(cfg_file, nameWolke, folderWolke, folderWolke, no_hist);     % 3rd input is where the hist files are, 4th input is where the output is saved
   
   % Keep only the particles, not the artifacts if any
   ids_prtcls = find(nube.pData.catPredict ~= 'Artifact');
   fieldnames = fields(nube.pData);
   for kf = 1:length(fieldnames)
      nube.pData.(fieldnames{kf}) = nube.pData.(fieldnames{kf})(ids_prtcls);
   end
   
   % Bring the dimensions to real positions in case that a spherical
   % illumination was used
   if nube.cfg.Ldivergent
      % Bring area, majsiz, minsiz, xpos, ypos and zpos back to pixels
      pixelsize   = nube.cfg.dx*nube.cfg.dy;
      nube.pData.area = nube.pData.area / pixelsize;
      nube.pData.majsiz = nube.pData.majsiz / sqrt(pixelsize);
      nube.pData.minsiz = nube.pData.minsiz / sqrt(pixelsize);
      nube.pData.xpos   = nube.pData.xpos / nube.cfg.dx + nube.parameter.Nx/2 + 1;
      nube.pData.ypos   = nube.pData.ypos / nube.cfg.dy + nube.parameter.Ny/2 + 1;
   
      % Get the magnification at every position
      magnifications = (nube.pData.zpos + nube.cfg.Ldivergent) / nube.cfg.Ldivergent;
      real_zpos = nube.pData.zpos ./ magnifications;   
   
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
end

% Plot histograms
subplot(2,2,1), histogram(real_zpos,20), title('Z position')
subplot(2,2,2), histogram(nube.pData.area,20), title('Area')
subplot(2,2,3), histogram(nube.pData.pDiam,20), title('Diameter (rescaled particle)')
subplot(2,2,4), histogram(nube.pData.asprat,20), title('Aspect ratio (majsiz/minsiz)')
pause

% Plot images showing the location of the particles
n_holos = max(nube.pData.holonum);
% Limits in x, y, z
fovx = (nube.parameter.Nx-1) * nube.parameter.dx / min(magnifications);
fovy = (nube.parameter.Ny-1) * nube.parameter.dy / min(magnifications);
z_range = [min(real_zpos), max(real_zpos)];
z_ints = uint8((real_zpos - z_range(1)) / (z_range(2) - z_range(1)) * 255);
% Initialize image
background = zeros(nube.parameter.Ny, nube.parameter.Nx, 3);
the_map = colormap('spring');
% Initialize spots
size_spots = 9;   % Use uneven number
[xx, yy] = meshgrid(linspace(-1,1,size_spots));
spot = sqrt(xx.^2 + yy.^2) <= 1;
[r_spot, c_spot] = find(spot);
% Take to relative positions
r_spot = r_spot - ceil((size_spots+1)/2);
c_spot = c_spot - ceil((size_spots+1)/2);
% Creating a movie
writerObj = VideoWriter('reconstructed_sequence.avi');
writerObj.FrameRate = 10;
open(writerObj);
figure
for kh = 1:n_holos
   id_prtcls = find(nube.pData.holonum == kh);
   if ~isempty(id_prtcls)
      values = ind2rgb(z_ints(id_prtcls),the_map);
      x_prtcl = nube.pData.xpos(id_prtcls);
      y_prtcl = nube.pData.ypos(id_prtcls);
      col_prtcl = round(x_prtcl *(nube.parameter.Nx - 1)/fovx  + (nube.parameter.Nx + 1)/2);
      row_prtcl = round(y_prtcl *(nube.parameter.Ny - 1)/fovy  + (nube.parameter.Ny + 1)/2);
      im = background;      
      for kp = 1:length(id_prtcls)
         % Get locations in image of the spot
         r_spot_prtcl = row_prtcl(kp) + r_spot;
         c_spot_prtcl = col_prtcl(kp) + c_spot;
         valids = r_spot_prtcl > 0 & r_spot_prtcl <= nube.parameter.Ny ...
            & c_spot_prtcl > 0 & c_spot_prtcl <= nube.parameter.Nx;
         r_spot_prtcl = r_spot_prtcl(valids);
         c_spot_prtcl = c_spot_prtcl(valids);
         i_spot_prtcl = sub2ind([nube.parameter.Ny, nube.parameter.Nx],r_spot_prtcl,c_spot_prtcl);
         % Apply colors to the image where the spots are
         for krgb = 1:3
            inds = i_spot_prtcl + (nube.parameter.Ny*nube.parameter.Nx) * (krgb-1);
            im(inds) = values(kp,1,krgb);
         end
      end
      % Show the image
      imshow(im)
      colormap('spring')
      title(sprintf('Particles in image %d/%d. Colobar limits : [%.2f,%.2f]mm',kh,n_holos,z_range(1)*1e3,z_range(2)*1e3))
      colorbar
      writeVideo(writerObj, getframe(gcf));
   end
end
close(writerObj)