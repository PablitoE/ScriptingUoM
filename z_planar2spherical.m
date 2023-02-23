% Get positions in z for direct propagation (planar wave) corresponding to
% a set of positions in z for a spherical wave set-up with 2 windows with
% different refractive index and thicknesses.
config_file = 'Config_Divergent.txt';

zs_plan = [24.65 29.5 42.1 69] *1e-3; % mm
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

zs_sph = zs_plan * L_eq ./ (L_eq+zs_plan);
fprintf('Equivalent zs in spherical propagation are: ')
fprintf('%.2f, ', zs_sph*1e3)
fprintf('mm\n')