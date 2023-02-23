% Script for performing a fine reconstruction around a set of particles.
% The location of maxima of tracers that are good for showing the
% presence of a bright spot are obtained for each particle and a prominency
% value of the peaks.
clear variables

% ClassificationData file with hand-labelled or predicted particles.
particle_detection = 'labelled'; % labelled or predicted
classdata_file = 'D:\HoloICE\Holograms\BrightSpot\coarse_reconstruction\oarse_reconstructio_cD.mat';
% Directory with holograms
holograms_dir = 'D:\HoloICE\Holograms\BrightSpot';
% Directory with hist files
hist_dir = 'D:\HoloICE\Holograms\BrightSpot\coarse_reconstruction';
% Directory for the new hist files. It will be in the same place as the
% hist directory of the coarse reconstruction
new_hist_dir = 'fine_reconstruction';
% Config file
config_file = 'D:\HoloICE\Holograms\BrightSpot\HoloBalloon_NYA_Wolke_coarse.cfg';
% Modifications to config file
mod_cfg = {'KzSampling', 0.5};
% Number of times the diameter of the particle that is reconstructed back and forwards in z.
n_diam = [3,3];

NumWorkers = 12;
reconHoloHist.refineSequenceOnParticles(config_file,classdata_file,particle_detection,new_hist_dir,mod_cfg,n_diam,1,NumWorkers);

% Get list of particles that might produce a bright spot
load(classdata_file,'temp1')
if strcmp(particle_detection,'labelled')
   indx_rnd_prtcl = temp1.label == 'Particle_round';
elseif strcmp(particle_detection,'predicted')
   indx_rnd_prtcl = temp1.predicted == 'Particle_round';
end
rnd_prtcls.ID = temp1.prtclID(indx_rnd_prtcl);
% Get location of particles and estimated diameters
parameters = {'equivdia','xpos','ypos','zpos'};
for p = parameters
   p = p{1};
   col_p = find(cellfun(@(x) strcmp(x,p),temp1.metricnames));
   rnd_prtcls.(p) = temp1.metricmat(indx_rnd_prtcl,col_p);
end
rnd_prtcls.n = length(rnd_prtcls.ID);
clear temp1
% Get list of hist files with fine reconstruction
path_fine_hist = fullfile(fileparts(hist_dir),new_hist_dir);
dir_out = dir(path_fine_hist);
files = {}; nf = 1;
for k=1:length(dir_out)
   if ~dir_out(k).isdir
      if endsWith(dir_out(k).name,'_hist.mat')
         files{nf} = fullfile(dir_out(k).folder,dir_out(k).name);
         nf = nf + 1;
      end
   end
end
% Get all the particles in the fine reconstruction
cfg = config(config_file);
cfg.current_sequence = 1;
ps = particles(cfg,files,numel(files),[],0);

n_prtcl_im = cellfun(@(x) numel(x), ps.holonum);   % Number of particles per image
fine_prtcls = struct([]);  % Initialization of list of fine reconstructed particles
filenames_in_ps = cell(length(ps.filenames),1);
% Remove path and _hist.mat from filenames in ps
for kf = 1:length(ps.filenames)
   [~,filenames_in_ps{kf}] = fileparts(ps.filenames{kf});
   filenames_in_ps{kf} = filenames_in_ps{kf}(1:end-5);
end
% Get the columns of the parameters of interest
for kp = 1:length(parameters) 
   p = parameters{kp};
   col_p(kp) = find(cellfun(@(x) strcmp(x,p),ps.prtclmetricvarnames));
end
% Get the columns of the traces of interest
traces_names = {'prtclvarlapltr','prtclprvarsobeltr'}; % bright spot detector and focus detector
for ktr = 1:length(traces_names) 
   tr = traces_names{ktr};
   col_tr(ktr) = find(cellfun(@(x) strcmp(x,tr),ps.prtcltracevarnames));
end
% Find the particles in the fine reconstruction that correspond to the
% coarsely detected particles.
for kp = 1:rnd_prtcls.n
   % Get the hologram filename that produced this particle. The index is the one in ps
   index_file = find(cellfun(@(x) startsWith(rnd_prtcls.ID{kp}{1},x),filenames_in_ps));
   prtcl_num_range = (sum(n_prtcl_im(1:index_file-1)) + 1) : (sum(n_prtcl_im(1:index_file))); % Range of particle numbers for this file
   % Find the closest particle in the fine reconstruction
   ratioDiam = ps.prtclmetrics(prtcl_num_range,col_p(1)) / rnd_prtcls.equivdia(kp);% < onePrtcl.eqsiz*0.2;
   distX = abs(rnd_prtcls.xpos(kp) - ps.prtclmetrics(prtcl_num_range,col_p(2)));%< min([onePrtcl.eqsiz/2000 0.1]);
   distY = abs(rnd_prtcls.ypos(kp) - ps.prtclmetrics(prtcl_num_range,col_p(3)));%< min([onePrtcl.eqsiz/2000 0.1]);
   distZ = abs(rnd_prtcls.zpos(kp) - ps.prtclmetrics(prtcl_num_range,col_p(4)));%< min([onePrtcl.eqsiz/50 1]);   
   dist = sqrt(distX.^2 + distY.^2 + distZ.^2);
   [~,idx_prtcl] = min(dist);
   idx_prtcl_total = prtcl_num_range(1) + idx_prtcl - 1;
   % Save traces and parameters
   for kpm = 1:length(parameters) 
      p = parameters{kpm};
      fine_prtcls(kp).(p) = ps.prtclmetrics(idx_prtcl_total,col_p(kpm));
   end
   for ktr = 1:length(traces_names)
      tr = traces_names{ktr};
      fine_prtcls(kp).(tr) = ps.prtcltraces{index_file}{idx_prtcl,col_tr(ktr)};
   end
   % Save ratio of diameters to check compatible particles
   fine_prtcls(kp).ratioDiam = ratioDiam(idx_prtcl);

   % Plot the image and the traces
   subplot(1,2,1), imagesc(abs(ps.prtclIm{index_file}{idx_prtcl})), colormap("bone")
   subplot(1,2,2), plot(fine_prtcls(kp).prtclprvarsobeltr,"DisplayName","Sobel-Focus"), hold on
   plot(fine_prtcls(kp).prtclvarlapltr,"DisplayName","Laplacian-BrightSpot")
   legend(), hold off
   pause()
end

