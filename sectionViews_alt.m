% Script for performing a fine reconstruction around a set of particles.
% The location of maxima of tracers that are good for showing the
% presence of a bright spot are obtained for each particle and a prominency
% value of the peaks.
clear variables

% ClassificationData file with hand-labelled or predicted particles.
particle_detection = 'labelled'; % labelled or predicted
classdata_file = 'D:\HoloICE\Holograms\BrightSpot\LN2\2023_01_31\sequence01\melting_freezing02_skip2\recon\recon_cD.mat';
% Auxiliar directory for the reconstruction
aux_dir = 'sectionViews_debug';
sectionViewsFilename = 'sV.mat';
update_info = false;  % Instead of recalculating the fine propagation or loading a previous result, just the structure with info is updated. Set as true or false.
% Config file
config_file = 'D:\HoloICE\Holograms\BrightSpot\LN2\2023_01_31\sequence01\melting_freezing02_skip2\recon\local_win.cfg';
% Sampling in z
KzSampling = 0.5;
% Number of times the diameter of the particle that is reconstructed back and forwards in x,y,z.
n_diam = [0.9,0.9;0.9,0.9;1,2.5];
NumWorkers = 12; % Number of processors for the reconstruction
n_random = [];  % Number of random particles to analyse. Use [] for using all particles.
max_z = []; % Maximum distance from hologram plane allowed for the particles. Use [] for using all particles.
min_num_z = 12; % Minimum amount of slices to save as section view of a particle. Use [] for using default.
normalise_traces = true;
toleranceMaxBrightSpot = 25;     % percentage of distance between bright spot and in-focus distances (EFL)
do_pause = false;
do_plot = false;
do_plot_transition = true;
do_scaleIms = true;
scaleImsRange = [0, 0.7];  % [0,1] would mean to use the min and max of all particles to scale all the images.
save_ims = false;
brightSpotBeforeAfter = 1; % -1 if bright spot is at a lower reconstruction distance (particle after the focal plane)
do_rescale = true;      % Rescale slices before calculating the new traces
bad_images_extract_from_plots = [];

% Plot slices and trace. Define sizes using normalized distances.
pos_trace_focus = [0.05, 0.575, 0.4, 0.15];
pos_trace_SL = [0.05, 0.8, 0.4, 0.15];
pos_sliceYZ = [0.01, 0.075, 0.48, 0.425];
pos_sliceXY = [0.5, 0.075, 0.48, 0.875];

OutputDir = fullfile(fileparts(classdata_file),aux_dir);
[status,~,~] = mkdir(OutputDir);
if ~status
   error("Directory output could not be created.")
end
file_output = fullfile(OutputDir,sectionViewsFilename);
if exist(file_output,'file')
   if update_info
      prtclSectionViews = reconHoloHist.getSectionViewsManyPrtcls(config_file,classdata_file,particle_detection,KzSampling,...
      n_diam,NumWorkers,'n_random',n_random,'max_z',max_z,'bright_spot',brightSpotBeforeAfter,'do_rescale',do_rescale,...
      'min_diam',1e-4,'update_info',file_output);
   save(file_output,"prtclSectionViews")
   else
      load(file_output,"prtclSectionViews")
   end
else
%    prtclSectionViews = reconHoloHist.getSectionViews(config_file,classdata_file,particle_detection,KzSampling,...
%       n_diam,NumWorkers,'n_random',n_random,'max_z',max_z,'bright_spot',brightSpotBeforeAfter,'do_rescale',do_rescale);
   prtclSectionViews = reconHoloHist.getSectionViewsManyPrtcls(config_file,classdata_file,particle_detection,KzSampling,...
      n_diam,NumWorkers,'n_random',n_random,'max_z',max_z,'bright_spot',brightSpotBeforeAfter,'do_rescale',do_rescale,...
      'min_diam',5e-5,'min_num_z',min_num_z,'keep_aspect_ratio',true);
   save(file_output,"prtclSectionViews")
end

% ANALYSIS
n_prtcls = size(prtclSectionViews, 1);
% Measures of lensing
laplacian_max = zeros(n_prtcls,1);
laplacian_max_bright_spot = zeros(n_prtcls,1);
laplacian_prominence = zeros(n_prtcls,1); % Calculated as max/min
laplacian_prominence_bright_spot = zeros(n_prtcls,1); % Calculated as maxBrightSpot/min
gini_max = zeros(n_prtcls,1);
gini_max_bright_spot = zeros(n_prtcls,1);
gini_prominence = zeros(n_prtcls,1); % Calculated as max/min
gini_prominence_bright_spot = zeros(n_prtcls,1); % Calculated as maxBrightSpot/min
ratio_max = zeros(n_prtcls,1);
ratio_max_bright_spot = zeros(n_prtcls,1);
ratio_prominence = zeros(n_prtcls,1); % Calculated as max/min
ratio_prominence_bright_spot = zeros(n_prtcls,1); % Calculated as maxBrightSpot/min

% Prepare the reconstruction of the in-focus particles after the section
% views became available. The search of the in-focus z position cannot be
% done while reconstructing if using multiprocess reconstruction.
propHandler = PropagationHandler(config_file);
for kp = 1:n_prtcls
   traceSobel = prtclSectionViews{kp,4}.traces.prtclprvarsobel;
   % Get the correct z position from the traces (they were calculated with
   % higher resolution and, probably, using the scaled patches).
   [~,zpos_index] = max(traceSobel);
   zs_trace = prtclSectionViews{kp,4}.zs;
   prtclSectionViews{kp,4}.xyz(3) = zpos_index;
   prtclSectionViews{kp,4}.zpos = zs_trace(zpos_index);
   % Update the in-focus image of the particles.
   propHandler.updateSlice(prtclSectionViews{kp,4}.hologram, prtclSectionViews{kp,4}.zpos);
   rows = prtclSectionViews{kp,4}.nys;
   cols = prtclSectionViews{kp,4}.nxs;
   prtclSectionViews{kp,3} = propHandler.slice(rows(1):rows(2),cols(1):cols(2));
end

% Prepare the scaling of the images
if do_scaleIms
   minmaxYZ = [min(cellfun(@(x) min(abs(x(:))),prtclSectionViews(:,2))),...
      max(cellfun(@(x) max(abs(x(:))),prtclSectionViews(:,2)))];
   scaleImYZ = minmaxYZ(1) + scaleImsRange * diff(minmaxYZ);
   minmaxXY = [min(cellfun(@(x) min(abs(x(:))),prtclSectionViews(:,3))),...
      max(cellfun(@(x) max(abs(x(:))),prtclSectionViews(:,3)))];
   scaleImXY = minmaxXY(1) + scaleImsRange * diff(minmaxXY);
end

for kp = 1:n_prtcls   
   % Get ranges of X, Y and Z
   nx = diff(prtclSectionViews{kp,4}.nxs) + 1;
   ny = diff(prtclSectionViews{kp,4}.nys) + 1;
   dx = prtclSectionViews{kp,4}.dxyz(1);
   dy = prtclSectionViews{kp,4}.dxyz(2);
   endX = dx * (nx-1);
   endY = dy * (ny-1);
   Zini_end = [prtclSectionViews{kp,4}.zs(1),prtclSectionViews{kp,4}.zs(end)];
   % The traces could come from the coarse reconstruction used for the
   % detection of round particles or from the finer reconstruction
   % performed for these section views (default).
   if isfield(prtclSectionViews{kp,4},'traces')
      zs_trace = prtclSectionViews{kp,4}.zs * 1e6; % um so that I can plot circles.
      traceSobel = prtclSectionViews{kp,4}.traces.prtclprvarsobel;
      traceLapla = prtclSectionViews{kp,4}.traces.prtclvarlapl;
      traceGini = prtclSectionViews{kp,4}.traces.prtclginilapl;
      traceRatio = prtclSectionViews{kp,4}.traces.prtclratiolapl;
      traceAmpg = prtclSectionViews{kp,4}.traces.prtclstdampg;
      traceCompg = prtclSectionViews{kp,4}.traces.prtclstdcompg;
   else
      % Get range of zs in the coarse reconstruction
      ind_coarse = find(prtclSectionViews{kp,4}.coarse.zs > Zini_end(1),1,"first"):find(prtclSectionViews{kp,4}.coarse.zs < Zini_end(2),1,"last");
      zs_trace = prtclSectionViews{kp,4}.coarse.zs(ind_coarse) * 1e6;
      % Prepare the traces. If lensing use trace for infocus (sobel) and
      % bright spot (laplacian). If not lensing use traces for detection of focus
      % (ampg and compg)
   %    if prtclSectionViews{kp,4}.islensing
      traceSobel = prtclSectionViews{kp,4}.coarse.prtclprvarsobeltr(ind_coarse);
      traceLapla = prtclSectionViews{kp,4}.coarse.prtclvarlapltr(ind_coarse);
      traceGini = prtclSectionViews{kp,4}.coarse.prtclginilapltr(ind_coarse);
      traceRatio = prtclSectionViews{kp,4}.coarse.prtclratiolapltr(ind_coarse);
      traceAmpg = prtclSectionViews{kp,4}.coarse.prtclstdampgtr(ind_coarse);
      traceCompg = prtclSectionViews{kp,4}.coarse.prtclstdcompgtr(ind_coarse);      
   end
   traceSobel_legend = "Var(Sobel of perimeter pixels)";
   traceLapla_legend = "Var(Laplacian of patch)";
   traceGini_legend = "Gini(abs(Lapl of thrIm))";
   traceRatio_legend = "MaxMid/Mean(AbsLapl)";
   traceAmpg_legend = "Std(Gradient of amplitude)";
   traceCompg_legend = "Std(Gradient of field)";   

   % Particle position to plot circle
   position = [(prtclSectionViews{kp,4}.xyz(1) - prtclSectionViews{kp,4}.nxs(1))*dx,...
      (prtclSectionViews{kp,4}.xyz(2) - prtclSectionViews{kp,4}.nys(1))*dy, prtclSectionViews{kp,4}.zpos] * 1e6; % um
   radius = prtclSectionViews{kp,4}.diam / 2 * 1e6;   % This value was taken from the rescaled particle if dynamic.addRescaleFeatures = true in config file.

   % Calculate measures of lensing with Laplacian
   laplacian_max(kp) = max(traceLapla);
   diam_um = prtclSectionViews{kp,4}.diam * 1e6;
   zBrightSpot = prtclSectionViews{kp,4}.zpos * 1e6 + brightSpotBeforeAfter * diam_um;
   zCloseToBrightSpot = (zs_trace > zBrightSpot - toleranceMaxBrightSpot/100*diam_um) & ...
      (zs_trace < zBrightSpot + toleranceMaxBrightSpot/100*diam_um);
   if ~any(zCloseToBrightSpot) % Case of very bad resolution in z or zBrightSpot out of range because of original bad zpos
      zCloseToBrightSpot(find(zs_trace < zBrightSpot,1,"last")) = true; % At least two z positions are considered as close to the bright spot
      zCloseToBrightSpot(find(zs_trace > zBrightSpot,1,"first")) = true;
   end
   laplacian_max_bright_spot(kp) = max(traceLapla(zCloseToBrightSpot));
   minLaplacian = min(traceLapla);
   laplacian_prominence(kp) = laplacian_max(kp) / minLaplacian;
   laplacian_prominence_bright_spot(kp) = laplacian_max_bright_spot(kp) / minLaplacian;
   % Calculate measures of lensing with Gini(Abs(Laplacian))
   gini_max(kp) = max(traceGini);
   gini_max_bright_spot(kp) = max(traceGini(zCloseToBrightSpot));
   minGini = min(traceGini);
   gini_prominence(kp) = gini_max(kp) / minGini;
   gini_prominence_bright_spot(kp) = gini_max_bright_spot(kp) / minGini;
   % Calculate measures of lensing with MaxMidPart/Mean(Abs(Laplacian))
   ratio_max(kp) = max(traceRatio);
   ratio_max_bright_spot(kp) = max(traceRatio(zCloseToBrightSpot));
   minRatio = min(traceRatio);
   ratio_prominence(kp) = ratio_max(kp) / minRatio;
   ratio_prominence_bright_spot(kp) = ratio_max_bright_spot(kp) / minRatio;

   if do_plot
      figure
      subplot('Position',pos_sliceYZ),
      %    imagesc(Zini_end*1e6,[0,endX]*1e6,abs(prtclSectionViews{kp,1}))
      %    xlabel('Z [um]'), ylabel('X [um]')
      %    axis equal
      %    viscircles([position(3), position(1)], radius);
      if do_scaleIms
         imYZ = imagesc(Zini_end*1e6,[0,endY]*1e6,abs(prtclSectionViews{kp,2}),scaleImYZ);
      else
         imYZ = imagesc(Zini_end*1e6,[0,endY]*1e6,abs(prtclSectionViews{kp,2}));
      end
      xlabel('Z [um]'), ylabel('Y [um]')
      axis equal
      viscircles([position(3), position(2)], radius,'LineStyle','--'); % Draw the particle
      xline(zBrightSpot - toleranceMaxBrightSpot/100*diam_um,'--g')
      xline(zBrightSpot + toleranceMaxBrightSpot/100*diam_um,'--g')
      zl = xlim;
      % Use limits of this YZ image to set the limits of the in focus image (XY)
      clim_global = imYZ.Parent.CLim;

      % Normalize traces to plot
      if normalise_traces
         normTrace = @(x) (x - min(x)) / (max(x) - min(x));
         traceSobel = normTrace(traceSobel);
         traceLapla = normTrace(traceLapla);
         traceGini = normTrace(traceGini);
         traceRatio = normTrace(traceRatio);
         traceAmpg = normTrace(traceAmpg);
         traceCompg = normTrace(traceCompg);
      end

      % Plot traces
      subplot('Position',pos_trace_focus), plot(zs_trace, traceAmpg, "DisplayName",traceAmpg_legend), hold on      
      plot(zs_trace, traceCompg, "DisplayName",traceCompg_legend)
      xlim(zl), legend('Location','southeast')
      xline(zBrightSpot - toleranceMaxBrightSpot/100*diam_um,'--g','HandleVisibility','off')
      xline(zBrightSpot + toleranceMaxBrightSpot/100*diam_um,'--g','HandleVisibility','off'), hold off      
      subplot('Position',pos_trace_SL), plot(zs_trace, traceSobel, "DisplayName",traceSobel_legend), hold on      
      plot(zs_trace, traceGini, "DisplayName",traceGini_legend),
      plot(zs_trace, traceLapla, "DisplayName",traceLapla_legend),
      plot(zs_trace, traceRatio, "DisplayName",traceRatio_legend),
      xlim(zl)
      legend('Location','southeast')
      xline(zBrightSpot - toleranceMaxBrightSpot/100*diam_um,'--g','HandleVisibility','off')
      xline(zBrightSpot + toleranceMaxBrightSpot/100*diam_um,'--g','HandleVisibility','off'), hold off

      % Plot in-focus image
      subplot('Position',pos_sliceXY)
      if do_scaleIms
         imagesc([0,endX]*1e6,[0,endY]*1e6,abs(prtclSectionViews{kp,3}),scaleImXY)
      else
         imagesc([0,endX]*1e6,[0,endY]*1e6,abs(prtclSectionViews{kp,3}),clim_global)
      end
      xlabel('X [um]'), ylabel('Y [um]')
      axis equal
      viscircles([position(1), position(2)], radius,'LineStyle','--');
      set(gcf,'Position',[680,300,1032,643]);
      colormap('bone')
      colorbar
      suptitle = sgtitle(sprintf("Particle from file: %s",prtclSectionViews{kp,4}.hologram));

      if save_ims
         figure_name_file = fullfile(OutputDir,sprintf("prtcl_%04d.svg",kp));
         saveas(gcf,figure_name_file)
         figure_name_file = fullfile(OutputDir,sprintf("prtcl_%04d.jpg",kp));
         saveas(gcf,figure_name_file)
      end
      if kp < n_prtcls
         if do_pause
            pause()
         end
         close(gcf)
      end
   end
end

if do_plot_transition
   nImages = 1:length(laplacian_max);
   nImages(bad_images_extract_from_plots) = [];   
   
   % Plot of possible lensing indicators. IF the analysed sequence is a
   % time sequence, it makes sense looking at them in a plot.
   figure
   subplot(1,2,1)
   plot(nImages,laplacian_max(nImages),'DisplayName','max(Laplacian)'), hold on
   plot(nImages,laplacian_max_bright_spot(nImages),'DisplayName','max(LaplacianBrightSpot)'), legend()
   title('Maximum Laplace value')
   subplot(1,2,2)
   plot(nImages,laplacian_prominence(nImages),'DisplayName','max/min(Laplacian)'), hold on
   plot(nImages,laplacian_prominence_bright_spot(nImages),'DisplayName','max/min(LaplacianBrightSpot)')
   legend()
   set(gcf,'Position',[385,550,1200,330])   

   figure
   subplot(1,2,1)
   plot(nImages,gini_max(nImages),'DisplayName','max(Gini)'), hold on
   plot(nImages,gini_max_bright_spot(nImages),'DisplayName','max(GiniBrightSpot)'), legend()
   title('Maximum Gini value (middle part of particles)')
   subplot(1,2,2)
   plot(nImages,gini_prominence(nImages),'DisplayName','max/min(Gini)'), hold on
   plot(nImages,gini_prominence_bright_spot(nImages),'DisplayName','max/min(GiniBrightSpot)')
   legend()
   set(gcf,'Position',[385,550,1200,330])  

   figure
   subplot(1,2,1)
   plot(nImages,ratio_max(nImages),'DisplayName','max(Ratio)'), hold on
   plot(nImages,ratio_max_bright_spot(nImages),'DisplayName','max(RatioBrightSpot)'), legend()
   title('Ratio of maximum Laplacian / mean Laplacian (middle part)')
   subplot(1,2,2)
   plot(nImages,ratio_prominence(nImages),'DisplayName','max/min(Ratio)'), hold on
   plot(nImages,ratio_prominence_bright_spot(nImages),'DisplayName','max/min(RatioBrightSpot)')
   legend()
   set(gcf,'Position',[385,550,1200,330])
  

   % Plot images of Laplacian at focus and bright spot locations
   % for k = 1:n_prtcls
   %    subplot(2,2,1), imagesc(abs(prtclSectionViews{k,3})),
   %    subplot(2,2,2), imagesc(abs(del2(abs(prtclSectionViews{k,3})))),colorbar
   %    subplot(2,2,3), imagesc(abs(prtclSectionViews{k,5})),
   %    subplot(2,2,4), imagesc(abs(del2(abs(prtclSectionViews{k,5})))),colorbar
   %    pause
   % end

   % Plot all the traces in a single figure
   figure
   particles = 1:15;
   traces_names = {'prtclprvarsobel','prtclvarlapl','prtclginilapl','prtclratiolapl','prtclstdampg'};
   n_traces = length(traces_names);
   np = length(particles);
   for kp = 1:np
      zs_trace = prtclSectionViews{kp,4}.zs;
      for kt = 1:n_traces
         subplot(1,n_traces,kt), plot(zs_trace,prtclSectionViews{kp,4}.traces.(traces_names{kt}))
         hold on
      end
   end
   for kt = 1:n_traces
      subplot(1,n_traces,kt), title(traces_names{kt})
   end
end