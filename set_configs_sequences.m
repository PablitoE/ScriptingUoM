% Script for creating config files from directories with sequences of
% images.
clear variables

% The directory with the sequences (folders). Two parts: root in this PC, changing is
% shared by JASMIN. Last directory in root is the name of the campaign
directory_sequences_root = 'D:\HoloICE\Holograms';
directory_sequences_changing = 'Pollen\11 - Pollen SVS\Ash_02';
% The directories with the sequences start with the same string
start_with = 'contrast_';
% Directories with the config files and results
directory_config_out_root = 'D:\HoloICE\MyHolosuite\Results';
directory_config_out_changing = 'Pollen\SVS\Ash_02';
camera = 'SVS-Vistek';  % Photron, Photron-pulsed, Optronis, SVS-Vistek, Simulation
L = 67.01e-3;
file_ext = 'png';
zmin = 28e-3;       % Range of z to sweep
zmax_guard = 4e-3;
Kzsampling = 3;
minPatchArea = 1;   % Minimum amount of pixels below threshold
numzs_ge = 8;
thresholdTuning = 0.5;  % Thresholding parameters
dilationMaskSize = 7;
useSAASM = false;
base_config = 'D:\HoloICE\MyHolosuite\Results\SimsSPIE\base_config.cfg';
nameFolderHist = 'multifocus';

jasmin_directory_sequences = fullfile_linux('/gws/nopw/j04/holography/pablo/holograms', directory_sequences_changing);
jasmin_config_out = fullfile_linux('~/Holosuite/Results/', directory_config_out_changing);

% Open the config file and have it ready to be saved again
cfg = config(base_config);
% [cfg_path, cfg_name, cfg_ext] = fileparts(base_config);   % Use the same name as the base config
% config_name = [cfg_name, cfg_ext];
config_name = 'local_win.cfg';
config_name_jasmin = 'jasmin.cfg';

% Update config
switch camera
   case 'Photron'
      dp = 20e-6;
      lambda = 532e-9;
   case 'Photron-pulsed'
      dp = 20e-6;
      lambda = 355e-9;   
   case 'Optronis'
      dp = 8e-6;
      lambda = 355e-9;
   case 'SVS-Vistek'
      dp = 2.2e-6;
      lambda = 355e-9;
   case 'Simulation'
      dp = 10e-6;
      lambda = 532e-9;
end
cfg.useSAASM = useSAASM;
cfg.dx = dp;
cfg.dy = dp;
cfg.lambda = lambda;
cfg.zMin = zmin;
cfg.zMax = L - zmax_guard;
cfg.KzSampling = Kzsampling;
cfg.Ldivergent = L;
cfg.hologram_filter = [cfg.hologram_filter(1:end-3), file_ext];
cfg.minPatchArea = minPatchArea;
% cfg.Padding = cfg.Padding(2:end-1);       % Remove ''  % SOLVED ELSEWHERE
the_rules = eval(cfg.dynamic.rules);
the_rules{cellfun(@(s)strcmp(s,'numzs'),the_rules(:,1)),3} = numzs_ge;
cfg.dynamic.rules = cellrules_to_evaluable(the_rules);
cfg.thresholdTuning = thresholdTuning;  % Thresholding parameters
cfg.dilationMaskSize = dilationMaskSize;
cfg.prefilters.amplitude{1,1}{1,2}.cutoffLenScale = dp;

% Parameters used in wolke
directory_config_out = fullfile(directory_config_out_root, directory_config_out_changing);
cfg.dynamic.NameInstrument = camera;
campaign = split(directory_sequences_root,'\');
campaign = campaign{end};
cfg.dynamic.NameCampaign = campaign;
pieces_sequences = split(cfg.sequences{1},',');
cfg.dynamic.predictor = fullfile(directory_config_out, [pieces_sequences{2}, '_tree.mat']);
cfg.dynamic.HistMinSize = 1e-6;
cfg.dynamic.IncludeTotal = 1; % True (use all kinds of particles-classes)
cfg.dynamic.RemoveRepeatedParticle = 1;   % Remove particles appearing at the same location
cfg.dynamic.SaveRepeatedParticles = 0;    % Save removed particles
cfg.dynamic.MaxCoincidenceScale = 0.00001;
% cfg.dynamic.zReScaling = 1;
% cfg.dynamic.maxScale = 1;
cfg.dynamic.loadCD = 0;
cfg.dynamic.earlyDiscardArtifacts = 1;
cfg.dynamic.rescaleParticles = 1;
cfg.dynamic.rescaleThresh = 0.45;
cfg.dynamic.iceLowestSize = dp/2;
% cfg.dynamic.waterMaxSize = dp*1000;

% Get the list of directories applying filter
directory_sequences = fullfile(directory_sequences_root, directory_sequences_changing);
files = dir(directory_sequences);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
% Get those that start with a particular string
flag_start_with_ok = zeros(length(subFolders),1, 'logical');
for k = 1:length(subFolders)
   flag_start_with_ok(k) = startsWith(subFolders(k).name, start_with,'IgnoreCase',true);
end
subFolders = subFolders(flag_start_with_ok);
Nseq = length(subFolders);
if Nseq == 0
   error("No contrast sequences were found.")
end
Nfiles_jobs = 0;

for kseq = 1:Nseq
   path = fullfile(subFolders(kseq).folder, subFolders(kseq).name);
   localTmp = fullfile(directory_config_out, subFolders(kseq).name);
   cfg.path = path;
   cfg.localTmp = localTmp;
   
   files = dir(path);
   extFlags = endsWith({files.name},file_ext);
   files = files(extFlags);
   if isempty(files)
      error('No images files found in sequence path. Try changing path or file_ext.')
   end
   Nfiles_jobs = Nfiles_jobs + length(files);
   T = struct2table(files); % convert the struct array to a table
	sortedT = sortrows(T, 'name'); % sort the table by 'DOB'
	files = table2struct(sortedT); % change it back to struct array if necessary
   sequences = [files(1).name, ':1:', files(end).name, ',', nameFolderHist];
   cfg.sequences = {sequences};
   
   sample_image = fullfile(files(1).folder, files(1).name);
   im = imread(sample_image);
   [Ny, Nx] = size(im);
   cfg.dynamic.Ny = num2str(Ny);
   cfg.dynamic.Nx = num2str(Nx);
   % cfg.dynamic.cD_data = fullfile(cfg.localTmp, [pieces_sequences{2}(1:end-1), '_cD.mat']);
   
   this_dir_output = fullfile(directory_config_out,subFolders(kseq).name);
   if ~exist(this_dir_output,'dir')
      mkdir(this_dir_output);
   end
   this_file_output = fullfile(this_dir_output,config_name);
   cfg.writefile(this_file_output);
   
   % JASMIN config file
   path_jasmin = fullfile_linux(jasmin_directory_sequences, subFolders(kseq).name);
   localTmp_jasmin = fullfile_linux(jasmin_config_out, subFolders(kseq).name);
   cfg.path = {path_jasmin, 'force'};
   cfg.localTmp = localTmp_jasmin;
   this_file_output = fullfile(this_dir_output,config_name_jasmin);
   cfg.writefile(this_file_output);
end
fprintf('Use up to %d jobs to process these sequences.\n', Nfiles_jobs)

function out = fullfile_linux(folder, file)   
   out = fullfile(folder, file);   
%    out = strrep(out,'\ ',' ');
   out(out == '\') = '/';
%    out = strrep(out,' ','\ ');
end

function s = cellrules_to_evaluable(C)
   s = '{';
   [n_rows,~] = size(C);
   for k = 1:n_rows
      if C{k,3} == int32(C{k,3})
         s = [s sprintf('''%s'',''%s'',%d',C{k,:})];
      else
         s = [s sprintf('''%s'',''%s'',%.1e',C{k,:})];
      end
      if k < n_rows
         s = [s ';'];
      end
   end
   s = [s '}'];
end