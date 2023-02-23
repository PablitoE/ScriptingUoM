% Get positions in z for direct propagation (planar wave) corresponding to
% a set of positions in z for a spherical wave set-up with 2 windows with
% different refractive index and thicknesses.
config_file = 'Config_Divergent.txt';

zs_sph = [20.4 29.4 42.7 50.8 55 57.3 63.7 64.7 66.8 67.5] *1e-3; % mm
str_input = get_inputs(config_file);
wavelength = str_input.wavelength;
L = str_input.L;
n1 = str_input.n1;
n2 = str_input.n2;
t1 = str_input.t1;
t2 = str_input.t2;

t1_eq = t1 / n1;
t2_eq = t2 / n2;
L_eq = L - t1 + t1_eq - t2 + t2_eq;

zs_plan = zs_sph * L_eq ./ (L_eq-zs_sph);
fprintf('Equivalent zs are: ')
fprintf('%.2f, ', zs_plan*1e3)
fprintf('mm\n')