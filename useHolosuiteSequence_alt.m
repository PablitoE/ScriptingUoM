open_carft = false;

% image_path = 'C:\Users\j05625pe\Documents\Work\Holograms\Spherical_Levitation\20210429\Seq2\Seq_end_12mm\image_drop.*.tiff';
% bgnd_path = 'C:\Users\j05625pe\Documents\Work\Holograms\Spherical_Levitation\20210429\Background\image_drop.*.tiff'; % [];
% PrepSeqHoloSuite(image_path,bgnd_path)
NumWorkers = 12;
if isempty(gcp("nocreate")) && NumWorkers > 1
   p = parpool('local', NumWorkers);
end
% Script for running Holosuite reconstruction of a sequence of holograms
if ispc
   if getComputerName() == 'l-u-as1002682'   % Waldemar PC
%       cfg_file = 'D:\HoloICE\MyHolosuite\Results\Glass beads\10 - New beads\Optronis 2um\Seq3\mode_0\contrast_seq0001\local_win.cfg';
%       cfg_file = 'D:\HoloICE\MyHolosuite\Results\Glass beads\10 - New beads\SVS Vistek 2um\Seq2\contrast_ims_scale0.36\local_win.cfg';
%       cfg_file = {'D:\HoloICE\MyHolosuite\Results\Glass beads\Photron\20210916_Microspheres2um\2um_inv1.5Ms_maxpower_C001H001S0001\contrast_seq0018\local_win.cfg';
%          'D:\HoloICE\MyHolosuite\Results\Glass beads\Photron\20210916_Microspheres2um\2um_inv1.5Ms_maxcarfpower_C001H001S0001\contrast_seq0022\local_win.cfg'};
%       cfg_file = 'D:\HoloICE\Holograms\Navitar 4k\2022_09_29\Upper\local_win.cfg'; % 'D:\HoloICE\Holograms\BrightSpot\LN2\round meets frosty_S0001\local_win.cfg'; % 
      cfg_file = 'D:\HoloICE\Holograms\BrightSpot\LN2\04_aggregating_timed\melting2_1110_to_1170_every_2\contrast_ims\local_win.cfg';
                  %'D:\HoloICE\Holograms\BrightSpot\LN2\04_aggregating_timed\melting1_410_to_439_every_2\contrast_ims\local_win.cfg';
                  %'D:\HoloICE\Holograms\BrightSpot\LN2\04_aggregating_timed\melting2_1110_to_1170_every_2\contrast_ims\local_win.cfg'};
      addpath('D:\HoloICE\MyHolosuite')
%       cfg_file = 'D:\HoloICE\MyHolosuite\Results\SimsSPIE\SAASM\mixed_size_num_particles_500\local_win_no_rules.cfg';
%       cfg_folder = 'D:\HoloICE\MyHolosuite\Results\SimsSPIE\SAASM';
%       if ~isempty(cfg_folder)
%          files = dir(cfg_folder);
%          dirFlags = [files.isdir];
%          dirFlags(1:2) = false;   % cases '.' and '..'
%          subFolders = files(dirFlags);
%          n_subfolders = length(subFolders);      
%          cfg_file = {};
%          cfg_file_name = 'local_win.cfg';
%          for kf = 1:n_subfolders
%             cfg_file{kf} = fullfile(subFolders(kf).folder,subFolders(kf).name,cfg_file_name);
%          end
%       end
   else
      cfg_file = 'C:\Work\Scripting\CompileHolosuite\Results\20210528\local_debug_win.cfg';
      addpath('C:\Work\Scripting\CompileHolosuite')
   end   
else
   cfg_file = '/home/pablo/Desktop/InHost/C_Work/Scripting/CompileHolosuite/Results/20210528/local_debug_ubu.cfg';
   addpath('/home/pablo/Desktop/InHost/C_Work/Scripting/CompileHolosuite')
end

% seqnum = 'image_drop.0._down1.png:1:image_drop.9._down1.png,whole';
% seqnum='C:\Users\j05625pe\Documents\Work\Holograms\Spherical_Levitation\20210429\Seq1';
% addpath('C:\Work\SmartGitHolosuite')

if iscell(cfg_file)
   num_seqs = length(cfg_file);
elseif ischar(cfg_file) || isstring(cfg_file)
   num_seqs = 1;
   cfg_file = {cfg_file};
end
order = 1:num_seqs;
         % [2, 1, 3:7]; % 100, 1000, 20, 200, 2000, 50, 500
for kseq = 1:num_seqs
   this_cfg_file = cfg_file{order(kseq)};
   reconHoloHist.reconSequence(this_cfg_file,1,NumWorkers);

   % Using seglen segnum, all images in the sequence are used for estimating
   % the background, good!
   % seglen = [8 8 8 8];
   % segnum = 2;
   % reconHoloHist.reconSequence(cfg_file,1,[],[],seglen,segnum);
   
   % index = '1';
   % n_job_array = '1';
   % jobID = '16';    % If less than 20, means amount of parallel Workers
   % profile -memory on
   % HolosuiteReconstructSequence(cfg_file,index,n_job_array,jobID)
   % profile viewer

   if open_carft
      [pathstr,name,ext] = fileparts(cfg_file);
      carft('config', fullfile(pathstr, 'whole', [name, ext]))
   end
end