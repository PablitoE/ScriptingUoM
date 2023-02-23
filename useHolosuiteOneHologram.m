% Script for running Holosuite reconstruction of a single hologram using
% the constant step-to-resolution ratio in the selection of zs position.
NumWorkers = 12;
p = gcp("nocreate");
if ~isempty(p) && p.NumWorkers ~= NumWorkers
   delete(p)
   p = gcp("nocreate");
end
if isempty(p) && NumWorkers > 1
   p = parpool('local',NumWorkers);
end
if ispc
   if getComputerName() == 'l-u-as1002682'   % Waldemar PC
%       cfg_file = 'D:\HoloICE\Spherical_Illumination_Holograms\base_config.cfg';
%       holo_file = 'simulated_image.png';
      cfg_file = 'D:\HoloICE\Holograms\BrightSpot\LN2\2023_01_31\sequence00\melting00_skip8\contrast_ims\local_win.cfg';
      holo_file = 'melting00_C001H001S0001000051.png';
%       cfg_file = 'D:\HoloICE\MyHolosuite\Results\Glass beads\Photron\20210916_Microspheres2um\2um_inv1.5Ms_maxpower_C001H001S0001\contrast_seq0022\local_win_single.cfg';
%       cfg_file = 'D:\HoloICE\MyHolosuite\Results\SimsSPIE\SAASM\mixed_size_num_particles_1000\local_win_less_rules.cfg';
%       holo_file = 'hologram_004.png'; % im_011 for seq0010
%       cfg_file = 'D:\HoloICE\MyHolosuite\Results\Splashing\local_win.cfg';
%       holo_file = 'cover_slip.png';
%       addpath('D:\HoloICE\MyHolosuite')
   else
      % Background image (using windows interfaces)
      cfg_file = 'C:\Work\Scripting\CompileHolosuite\Results\20210329\Background\local_analysis_win.cfg';
      holo_file = 'image_drop._down1_aver.png';

      % Image with a single particle
      cfg_file = 'C:\Work\Scripting\CompileHolosuite\Results\20210429\Seq1\local_debug_win.cfg';
      holo_file = 'image_drop.0.png';
   end
else
   % Image with a single particle
   cfg_file = '/home/pablo/Desktop/InHost/C_Work/Scripting/CompileHolosuite/Results/20210429/Seq1/local_debug_ubu.cfg';
   holo_file = 'image_drop.0.png';
end
% addpath('C:\Work\SmartGitHolosuite')
% zs = sample_z_sph(10e-3,80e-3,500,82e-3,3e-2,355e-9);

if iscell(holo_file)
   for kim = 1:length(holo_file)
      hf = holo_file{kim};
      times = reconHoloHist.reconHoloHistAndSave(cfg_file, hf);
   end
else
   times = reconHoloHist.reconHoloHistAndSave(cfg_file, holo_file);
end

carft('config',cfg_file)