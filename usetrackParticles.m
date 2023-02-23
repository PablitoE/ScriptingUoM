% Script for learning how to use trackParticles class
clear variables

results_directory = 'D:\HoloICE\MyHolosuite\Results\Glass beads\10 - New beads\Optronis 2um\Seq3\contrast_seq0001';
files_in_directory = dir(results_directory);

% Get config file and classification file
k = 1;
cfg_file_found = false;
class_file_found = false;
while k <= length(files_in_directory) && ~(cfg_file_found  && class_file_found)
   current_file = files_in_directory(k);
   if endsWith(current_file.name, '.cfg') && strcmp(current_file.name(1:5), 'local')
      cfg_file = fullfile(current_file.folder, current_file.name);
      cfg_file_found = true;
   end
   if endsWith(current_file.name, 'cD.mat')
      class_file = fullfile(current_file.folder, current_file.name);
      class_file_found = true;
   end
   k = k + 1;
end

% carft('config',cfg_file)
% load(class_file)

% % Get a hist file
% for k = 1:length(files_in_directory)
%    current_file = files_in_directory(k);
%    if length(current_file.name) > 4 && strcmp(current_file.name(end-7:end), 'hist.mat')
%       load_struct = load(fullfile(current_file.folder, current_file.name));
%       pStats = load_struct.pd.getpStats();
%    end
% end

% pStats is a structure with all the pStats from the files concatenated.
% zs, Nx, Ny are the only fields that are not concatenated as all files
% have the same information.
% The field holonum indicates from which hologram it came.
% The constructor of particles or addfiles could be used.

tracker = trackParticles(cfg_file, class_file);
tracker.makeTracks();