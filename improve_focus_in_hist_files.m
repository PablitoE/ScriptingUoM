% Script for improving the detected focus of particles. A particledata
% object is created and saved in new hist files. Also a new cD file is
% saved which has the information of the updated round particles.
clear variables

particle_detection = 'labelled'; % labelled or predicted
classdata_file = 'D:\HoloICE\ICE-D_mainz\holograms\warm\seq12\seq12_cD.mat'; % ClassificationData file with hand-labelled or predicted particles.
config_file = 'D:\HoloICE\ICE-D_mainz\holograms\warm\seq12\kotelett.cfg';
NWorkers = 1;
MinPeakDistance = 4;    % Parameter for looking for peaks in the traces
nKeepPatches = 120;     % Amount of patches to keep from every particle
output_dir = 'corrected_hists';
ValidatedByHandFilename = 'validatedFocuses';
saveOnlyRoundParticles = true;
do_plot = 0;      % 0: no plot, 1: just refocused and pause, 2: all, 3: all and pause in refocused
rescaleThreshold = 0.6; % Rescaled threshold parameter for getting scores of possible in-focus patches.
validatingByHand = false; % Longer times to show plots so that we save records of properly validated focuses.
JustUseChosenFromValidated = true; % Recalculate cD and hist files but get the chosen fragment from the ValidatedPrtcls file.

% Define inline functions
norm1d = @(x) (x - min(x)) ./ (max(x) - min(x));
std_med = @(x) mean((x(:) - median(x(:))).^2);

% Load the config file
if ischar(config_file)
   cfg = config(config_file);
   fprintf('Config file opened: %s\n', config_file);
end
[holofiles, savedir] = cfg.seq_list;

% Choose the number of workers
p = gcp('nocreate'); % Get the pool size if there is one
if isempty(p), poolsize = 0; else, poolsize = p.NumWorkers; end
% Setup the parpool as needed
if NWorkers <= 1 && poolsize, delete(gcp('nocreate')); end % Close the pool if need be
if ~isfinite(NWorkers) && poolsize == 0    % Open the default pool if size is not specified
   parpool('local');
elseif NWorkers > 1 && poolsize ~= NWorkers  % Reopen pool to desired size
   p = gcp('nocreate'); % Shut the parpool down if there is one
   if ~isempty(p), delete(gcp('nocreate')); end
   parpool('local',NWorkers);% Start the parpool with the desired number of workers
end

% Set the waiting times
waitNotRefocusing = 0.5;
validatingByHand = validatingByHand && do_plot == 3;   % Need to be watching carefully to be really validating decisions by hand.
if validatingByHand, waitNotRefocusing = 1; end

% Information about the setup
Nx = str2double(cfg.dynamic.Nx);
Ny = str2double(cfg.dynamic.Ny);
D = min(Nx * cfg.dx, Ny * cfg.dy);  % Minimum dimension of (magnified) sensor
% Get list of reconstructed particles
hash_cD = DataHash(classdata_file,'MD5','hex','file');
cD = load(classdata_file,'temp1');
cD = cD.temp1;
if ~issorted(cD.holonum)
   error("The cD file should be ordered in holonums.")
end
if strcmp(particle_detection,'labelled')
   indx_rnd_prtcl = cD.label == 'Particle_round';
elseif strcmp(particle_detection,'predicted')
   indx_rnd_prtcl = cD.predicted == 'Particle_round';
end
idx_artifacts = find(cD.label == "Artifact");
n_metrics_cD = length(cD.metricnames);
% If the input n_random was given a value, just work with some particles
kps_in_cD = find(indx_rnd_prtcl);
n_input_particles = length(kps_in_cD);
hist_files = cD.filenames;
% Directory with cData file has also the hist files
[dir_data, cD_filename, cD_ext] = fileparts(classdata_file);
n_holofiles = length(holofiles);
length_ext = length(cfg.hologram_filter) - 1;

% Prepare the output directory
output_path = fullfile(dir_data,output_dir);
if ~exist(output_path,"dir")
   mkdir(output_path);
end
savingcD_file = fullfile(output_path,[cD_filename, cD_ext]);

% Prepare names of hologram files
holofilenames = holofiles;
for kholo = 1:n_holofiles
   [~,holofilenames{kholo},~] = fileparts(holofilenames{kholo});
end

string = sprintf('Prepared %d particles from cD file %s.\n', ...
   n_input_particles,classdata_file);

% Prepare reading of sobel data from patches of interest to get score
IdxVarSobel = find(strcmp('prtclprvarsobel',reconHoloHist.prtclvarnames));

% Initializations
MsgError = [];
current_hist_file_index = 0;
current_holonum_loaded = 0;
changes_made_to_pd = false;
nPrtclsIncD = length(cD.prtclID);
remainigPrtclsIncD = true(nPrtclsIncD,1); % Keep track of the particles that haven't been analysed yet
namesOfFocusingScores_ = {'scores_mean','ratioscores','scores','difscores','scores_rescaled',...
   'scores_mean_rescaled','difscores_rescaled','ratioscores_rescaled','prtclPerimeterSobel'};
NfocusingScores = length(namesOfFocusingScores_);

% Read validated particles in savingcD_file
if ~endsWith(ValidatedByHandFilename,'.mat'), ValidatedByHandFilename = [ValidatedByHandFilename '.mat']; end
savingValidated = fullfile(output_path,ValidatedByHandFilename);
if exist(savingValidated,"file")
   load(savingValidated,"validatedPrtcls","focusingScores","chosenFocus","namesOfFocusingScores")
   if ~all(strcmp(namesOfFocusingScores,namesOfFocusingScores_))
      error("New focusing scores are available. The ValidatedByHand file should be removed before rerunning.")
   end
else
   validatedPrtcls = false(nPrtclsIncD,1);
   focusingScores = cell(n_input_particles,NfocusingScores);
   chosenFocus = zeros(n_input_particles,1);
   namesOfFocusingScores = namesOfFocusingScores_;
end
clear namesOfFocusingScores_;

if do_plot == 0
   h = waitbar(0,"Correcting focus...");
end
try
for kp = 1:n_input_particles
   kpcD = kps_in_cD(kp);
   % Continue if this particle's focus point has been validated and saved
   if ~validatedPrtcls(kpcD) || JustUseChosenFromValidated

      hist_file = hist_files{cD.holonum(kpcD)}; % withoutpath
      index_file = find(cellfun(@(x) startsWith(hist_file,x(1:end-length_ext)),holofiles));
      % The information about the focusing is saved in the hist file
      if current_holonum_loaded ~= cD.holonum(kpcD)
         % Only load it if the previously loaded hist file is different
         current_holonum_loaded = cD.holonum(kpcD);
         hist_file = fullfile(dir_data,hist_file);
         pd = load(hist_file);
         if isfield(pd,'pStats')
            pd = particledata(pd.pStats,pd.b,pd.times,[],holofilenames{index_file});
         else
            pd = pd.pd;
         end
         n_metrics = length(pd.metricnames);
         n_objects = pd.Nobjects;
         % Get the particles' number in this hist file that are present
         % in the classified particles (which are not all necessarily
         % round particles)
         [~,ia,ib] = cD.intersect(pd,remainigPrtclsIncD);
         artifactsInPD = ib(ismember(ia,idx_artifacts));
         roundParticlesInPD = ib(ismember(ia,kps_in_cD));
         remainigPrtclsIncD(ia) = false;
      end
      % kp points to one round particle among the classified particles,
      % which can now be linked to a detected particle in the hist file.
      kh = ib(ia == kpcD);

      % Work with the traces prtclstdampgtr and prtclstdcompgtr
      prtclstdampgtr = pd.gettrace("prtclstdampgtr",kh);
      prtclstdcompgtr = pd.gettrace("prtclstdcompgtr",kh);
      % zLIndampg, zLIndcompg : a priori positions of best peak
      % locs_amp,locs_comp,start_fragments,end_fragments: description of all possible peaks and fragments where the particle could live in
      % cnt_fragment: fragment with more chances of being the best
      % n_fragments: amount of possible fragments
      [zLIndampg,zLIndcompg,locs_amp,locs_comp,start_fragments,end_fragments,cnt_fragment,n_fragments] = reconHoloHist.getFragmentsWithFoci(prtclstdampgtr,prtclstdcompgtr,MinPeakDistance);      
      if n_fragments > 1
         % Keep only one fragment with the best score of the patch (more pixels well below threshold)
         locs_PatchOfInterest = round((locs_comp + locs_amp)/2);
         % Retrieve patches
         PatchOfInterest = pd.getprtclImtrObjAt(kh,locs_PatchOfInterest);
         % Get the list of used thresholds. Normally there is a single
         % threshold for all the objects coming from a hologram but this could
         % change in the future and the information is available
         threshs = pd.gettrace("prtclthreshtr",kh);
         threshs = num2cell(threshs(locs_PatchOfInterest));
         % Get a score for each of the images
         scores_mean = cellfun(@(x,thr) reconHoloHist.patch_score(~isnan(x),abs(x),thr,@mean),PatchOfInterest,threshs);
         scores_std = cellfun(@(x,thr) reconHoloHist.patch_score(~isnan(x),abs(x),thr,std_med),PatchOfInterest,threshs);
         scores = cellfun(@(x,thr) reconHoloHist.patch_score(~isnan(x),abs(x),thr,@median),PatchOfInterest,threshs);
         antiscores = cellfun(@(x,thr) reconHoloHist.patch_score(~isnan(x),thr - abs(x),thr,@median),PatchOfInterest,threshs);
         difscores = scores - antiscores;
         %       ratioscores = scores ./ antiscores;
         ratioscores = scores_mean ./ scores_std;
         nPOI = length(PatchOfInterest);
         scores_rescaled = zeros(nPOI,1); antiscores_rescaled = scores_rescaled;
         scores_mean_rescaled = scores_rescaled; scores_std_rescaled = scores_rescaled;
         prtclPerimeterSobel = zeros(nPOI,1);
         %       figure
         for kPOI = 1:nPOI
            [thrslice,slice] = particles.rescaleImage(PatchOfInterest{kPOI},rescaleThreshold,cfg.closeGapsRad,cfg.shouldFillHoles);
            %          [r,c] = find(thrslice);
            %          imagesc(abs(slice)), hold on, plot(r,c,'xr'), hold off
            scores_rescaled(kPOI) = reconHoloHist.patch_score(~isnan(slice),abs(slice),rescaleThreshold,@median);
            antiscores_rescaled(kPOI) = reconHoloHist.patch_score(~isnan(slice),rescaleThreshold - abs(slice),rescaleThreshold,@median);
            scores_mean_rescaled(kPOI) = reconHoloHist.patch_score(~isnan(slice),rescaleThreshold - abs(slice),rescaleThreshold,@mean);
            scores_std_rescaled(kPOI) = reconHoloHist.patch_score(~isnan(slice),rescaleThreshold - abs(slice),rescaleThreshold,std_med);
            % Get the traces values for these rescaled patches
            PixelIdxList = find(~isnan(slice)); slice_size = size(slice);
            prtcl_data = reconHoloHist.getprtcldaten(PixelIdxList, slice(PixelIdxList), rescaleThreshold, cfg.dilationMaskType,...
               cfg.dilationMaskSize, cfg.closeGapsRad, cfg.shouldFillHoles, slice_size(2), slice_size(1));
            prtclPerimeterSobel(kPOI) = prtcl_data(IdxVarSobel);
         end
         difscores_rescaled = scores_rescaled - antiscores_rescaled;
         %       ratioscores_rescaled = scores_rescaled ./ antiscores_rescaled;
         ratioscores_rescaled = scores_mean_rescaled ./ scores_std_rescaled;

         % Save the scores
         for ks = 1:NfocusingScores
            focusingScores{kp,ks} = eval(namesOfFocusingScores{ks});
         end

         if ~(validatedPrtcls(kpcD) && JustUseChosenFromValidated && chosenFocus(kp)) % Only update the chosenFocus if it's required
            [~,chosenFocus(kp)] = max(scores); % It looks like the best thing is to use the direct score with the median.
         end
         zLInd = locs_PatchOfInterest(chosenFocus(kp));
         zLIndampg = locs_amp(chosenFocus(kp));
         zLIndcompg = locs_comp(chosenFocus(kp));
      else
         zLInd = round(mean([zLIndampg,zLIndcompg]));
      end
      new_focus = pd.zLInd(kh) ~= zLInd;

      % Plotting
      if ismember(do_plot,[2,3]) || (do_plot == 1 && new_focus)
         subplot(1,3,1)
         hold off, plot(prtclstdcompgtr), hold on
         plot(prtclstdampgtr), title(sprintf("Particle %d of %d. Red vertical line is the previous focus position.",kp,n_input_particles))
         xline(pd.zLInd(kh),'r')
         xline(zLInd,'g')
%          plots_score = norm1d([scores, scores_mean, scores_std, difscores, scores_rescaled, scores_mean_rescaled, scores_std_rescaled, difscores_rescaled]);
%          subplot(1,4,2), plot(locs_PatchOfInterest, plots_score')
%          legend('scores', 'scores_mean', 'scores_std', 'difscores', 'scores_rescaled', 'scores_mean_rescaled', 'scores_std_rescaled', 'difscores_rescaled')
         subplot(1,3,2), imagesc(abs(cD.prtclIm{kpcD})), title("Old in-focus patch")
      end

      % New focus point case
      if new_focus
         % Update this object in pd with a new focus point.
         pd.zLInd(kh) = zLInd;
         % Update metrics from traces
         metric_done = zeros(n_metrics,1,'logical');  % Keep track of metrics that were updated
         metricsUnchanged = {'numzs','phfl'};   % Metrics that don't change with a change of the focus position
         for km = 1:n_metrics
            if ismember(pd.metricnames{km},metricsUnchanged)   % If nothing to do, continue
               metric_done(km) = true;
               continue
            end
            % The trace corresponding to the metric can take one of two names:
            idx_tracetr = find(strcmp([pd.metricnames{km} 'tr'],pd.tracenames));
            idx_prtcltracetr = find(strcmp(['prtcl' pd.metricnames{km} 'tr'],pd.tracenames));
            if isempty(idx_prtcltracetr) && isempty(idx_tracetr)
               idx_trace = [];
            elseif isempty(idx_prtcltracetr) && ~isempty(idx_tracetr)
               idx_trace = idx_tracetr; tracename = [pd.metricnames{km} 'tr'];
            elseif ~isempty(idx_prtcltracetr) && isempty(idx_tracetr)
               idx_trace = idx_prtcltracetr; tracename = ['prtcl' pd.metricnames{km} 'tr'];
            else
               error("I can't resolve having two traces corresponding to a single metric.")
            end
            if ~isempty(idx_trace)  % Update the metric if the trace was found
               the_trace = pd.gettrace(tracename,kh);
               pd.metricmat(kh,km) = the_trace(zLInd);
               metric_done(km) = true;
            end
         end
         % Update the other metrics
         eqsiz = mean([pd.getmetric('majsiz',kh), pd.getmetric('minsiz',kh)]);
         particular_metrics = {'zLIndampg', zLIndampg; 'zLIndcompg', zLIndcompg; 'zLInd', zLInd; 'eqsiz', eqsiz};
         for kpm = 1:size(particular_metrics,1)
            idx_pm = find(strcmp(pd.metricnames,particular_metrics{kpm,1}));
            pd.metricmat(kh,idx_pm) = particular_metrics{kpm,2};
            metric_done(idx_pm) = true;
         end
         unresolved_metrics = pd.metricnames(~metric_done);
         if ~isempty(unresolved_metrics)
            error("Uploading hist file. Couldn't resolve metrics: %s", strjoin(unresolved_metrics))
         end
         % Update prtclID
         pd.prtclID{kh} = reconHoloHist.makeParticleID(holofilenames{index_file}, ...
            pd.getmetric('xpos',kh), pd.getmetric('ypos',kh), ...
            pd.getmetric('zpos',kh), pd.getmetric('eqsiz',kh));

         % Update the cD too. Changes in prtclID, metricmat and prtclIm
         cD.prtclID{kpcD} = pd.prtclID{kh};
         cD.prtclIm{kpcD} = pd.getprtclImtrObjAt(kh,zLInd);
         cD.prtclIm{kpcD} = cD.prtclIm{kpcD}{1};
         metric_done = false(n_metrics_cD,1);
         for kmcd = 1:n_metrics_cD  % Start with metrics available in pd
            idx_metricInpd = find(strcmp(cD.metricnames{kmcd},pd.metricnames));
            if ~isempty(idx_metricInpd)
               cD.metricmat(kpcD,kmcd) = pd.metricmat(kh,idx_metricInpd);
               metric_done(kmcd) = true;
            elseif startsWith(cD.metricnames{kmcd},'tph')
               metric_done(kmcd) = true; % The Peak Height metrics don't change with a change in the focus position
            end
         end
         % Continue updating metrics not included in pd
         equivdia = (4/pi*pd.getmetric('area',kh)).^(1/2);
         fakepStats = struct; % Use a fake pStats to get other
         particular_metrics = {'equivdia', equivdia};
         for kpm = 1:size(particular_metrics,1)
            idx_pm = find(strcmp(cD.metricnames,particular_metrics{kpm,1}));
            cD.metricmat(kpcD,idx_pm) = particular_metrics{kpm,2};
            metric_done(idx_pm) = true;
         end
         unresolved_metrics = cD.metricnames(~metric_done);
         % Leave other metrics as they are (see reconHoloHist.addMoreMetrics)
         %       if ~isempty(unresolved_metrics)
         %          error("Uploading cD file. Couldn't resolve metrics: %s", strjoin(unresolved_metrics))
         %       end
         if do_plot
            subplot(1,3,3), imagesc(abs(cD.prtclIm{kpcD})), title("New in-focus patch");
            if do_plot == 2, pause(1), else, pause(), end
         end
      elseif ismember(do_plot,[2,3])
         subplot(1,3,3), imagesc([]), title("New in-focus patch")
         pause(waitNotRefocusing)
      end
      % Reduce the size of the traces by keeping only the neighbors of the
      % focus. This will change the properties: tracemat, numzsind, numxsind,
      % numysind, prtclimvec, prtclthreshimvec; and the metrics numzs and phfl
      [pd,changes_made] = pd.reduceTracesOfObj(kh,nKeepPatches);
      if changes_made % Update metrics in cD that changed due to the trace trimming
         particular_metrics = {'numzs','zLIndampg','zLIndcompg','zLInd','phfl'};
         for kpm = 1:size(particular_metrics,1)
            idx_pm_pd = find(strcmp(pd.metricnames,particular_metrics{kpm}));
            idx_pm_cD = find(strcmp(cD.metricnames,particular_metrics{kpm}));
            cD.metricmat(kpcD,idx_pm_cD) = pd.metricmat(kh,idx_pm_pd);
         end
         % The Peak Height metrics change with a trimmed trace
         for kmcd = 1:n_metrics_cD
            if startsWith(cD.metricnames{kmcd},'tph')
               tracename = ['prtcl' cD.metricnames{kmcd}(4:end) 'tr'];
               if sum(strcmp(tracename,pd.tracenames)) == 1
                  cD.metricmat(kpcD,kmcd) = particles.estPeakHeight(pd.gettrace(tracename,kh));
               else
                  error("There is no trace information to update %s of cD.",cD.metricnames{kmcd})
               end
            end
         end
      end
      changes_made_to_pd = changes_made || changes_made_to_pd;
      if validatingByHand
         validatedPrtcls(kpcD) = true;
      end
   end % Only done if not a validated particle

   % Save new hist file if finished with these (and this holo hist has been opened)
   if current_holonum_loaded > 0 && (kp == n_input_particles || (current_holonum_loaded == cD.holonum(kpcD) && ~ismember(kps_in_cD(kp+1),ia)))
      if saveOnlyRoundParticles
         changes_made_to_pd = true;
         otherThanRoundParticles = setdiff(1:n_objects,roundParticlesInPD);
         pd = pd.removeObjects(otherThanRoundParticles);
      end
      if changes_made_to_pd
         savinghist_file = fullfile(output_path,hist_files{current_holonum_loaded});
         pd.prtclthreshimvec = []; % Save space as this can be reconstructed easily.
         save(savinghist_file,'pd','-v7.3');
      end
   end
   if do_plot == 0
      waitbar(kp/n_input_particles,h)
   end   
end
close(h)
catch MsgError
   if validatingByHand      
      lastValidatedKp = input("Last validated particle: ");
      kp_ignore = lastValidatedKp+1 : kp;
      for kpi = 1:length(kp_ignore), for ks = 1:NfocusingScores, focusingScores{kp_ignore(kpi),ks} = []; end, end
      chosenFocus(kp_ignore) = 0; kpcd_ignore = kps_in_cD(kp_ignore);
      validatedPrtcls(kpcd_ignore) = false;
      save(savingValidated,"validatedPrtcls","focusingScores","chosenFocus","namesOfFocusingScores")
   end
end
temp1 = cD;
clear cD
save(savingcD_file,"temp1")
clear temp1
if validatingByHand && kpcD == kps_in_cD(n_input_particles) && validatedPrtcls(kpcD)
   save(savingValidated,"validatedPrtcls","focusingScores","chosenFocus","namesOfFocusingScores")
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = gcp('nocreate'); % Shut the parpool down if there is one
if ~isempty(p), delete(gcp('nocreate')); end

if ~isempty(MsgError)
   disp("Error catched and managed.")
   rethrow(MsgError)
end