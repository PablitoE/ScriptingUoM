% Save config file for using with sweep_pars_SIFT_JASMIN
clear variables
z_file = 'C:\Work\Scripting\results\Calibration USAF\weigthed_average\z_positions.mat';
ims_folder = 'C:\Users\j05625pe\Documents\Work\Holograms\Glass beads\Photron - 2um\20210812\Calibration USAF';
output_file = '.\results\output_debug';
dx = 20e-6;
wavelength = 532e-9;
Np = 20; Ne = 15; Nt = 10;
Nh = 5;
bounds_p = [3, 50];
bounds_e = [3, 30];
bounds_t = [0.9, 4];
bounds_h = [0.02, 0.2];

save('config_SIFT_sweep.mat')