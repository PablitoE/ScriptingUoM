clear variables
close all
mm = 1e-3;
um = 1e-6;
nm = 1e-9;

L_guess = 90*mm;  % Guessed distance from sensor to source
% zs = [21.35, 29.50, 65.05, 67.90]*mm;  % Distances from sensor when using a L guessed (thicknesses 2.3 mm 11.8mm)
% zs = [20.14, 29.26, 57.1, 63.55]*mm; % 20201117 Distances from sensor when using a L guessed of 80e-3 (thicknesses 5 mm 12mm)
zs = [11, 21.89, 28.19, 37.93]*mm;

if length(zs) == 4
   n_windows_estimation = 2;
elseif length(zs) == 2
   n_windows_estimation = 2;
else
   error('Amount of zs for estimation not supported.')
end
      
nFS = 1.4761; % Fused Silica at 355 nm
n1 = nFS;
t1_data = 5 *mm;  % Thickness of window close to the source used for getting data, the actual thickness was 5 mm
t1 = 5*mm;
n2 = nFS;
t2_data = 12 *mm; % Thickness of window close to the sensor
t2 = 12*mm;
diag_sensor = 35.48 *mm;   % SVS-Vistek hr120mcx
dx_sensor = 2.2 *um;       % SVS-Vistek hr120mcx
wavel = 355 *nm;

t1_eq_data = t1_data / n1;  % Effective thicknesses
t1_eq = t1 / n1;
t2_eq_data = t2_data / n2;
t2_eq = t2 / n2;
L = L_guess - t1_data + t1_eq_data - t2_data + t2_eq_data;    % Effective L for getting data (zs)
fprintf('Reconstruction effective L = %.2f mm\n',L/mm)

z_recon = zs * L ./ (L-zs);   % Distances of reconstruction that were actually used for obtaining the values of z for the windows' interfaces

% Getting the L that minimices the error of zs
if n_windows_estimation == 2
   syms z1 z2 z3 z4 L t1eff t2eff
   Q1 = (z2*L/(z2+L) - z1*L/(z1+L) - t2eff)/t2eff; 
   Q2 = (z4*L/(z4+L) - z3*L/(z3+L) - t1eff)/t1eff;
   Q = Q1^2 + Q2^2;
   Q_eval = subs(Q,[z1 z2 z3 z4 t1eff t2eff],[z_recon t1_eq t2_eq]);
else
   syms z1 z2 L t1eff
   Q1 = (z2*L/(z2+L) - z1*L/(z1+L) - t2eff)/t2eff; 
   Q = Q1^2;
   Q_eval = subs(Q,[z1 z2 t2eff],[z_recon t2_eq]);
end
dQdL = diff(Q_eval,L);
L_eff_solved = solve(dQdL==0,L);
L_eff = double(vpa(L_eff_solved));
L_eff = L_eff(imag(L_eff) == 0);    % Remove complex solutions
L_eff = L_eff(L_eff > 0);    % Remove negative or zero solutions
fprintf('Effective L = %.2f mm.\n',L_eff/mm)

Total_L = L_eff + t1 - t1_eq + t2 - t2_eq;
Q1err = double(subs(Q1,[z1 z2 t2eff L],[z_recon(1:2) t2_eq L_eff])) * t2_eq;
fprintf('Real L = %.2f mm.\n',Total_L/mm)
fprintf('Error in thickness of window 1 = %.2f mm.\n',Q1err/mm)
if n_windows_estimation == 2
   Q2err = double(subs(Q2,[z3 z4 t1eff L],[z_recon(3:4) t1_eq L_eff])) * t1_eq;
   fprintf('Error in thickness of window 2 = %.2f mm.\n',Q2err/mm)
end

% Distances sensor to windows (first interface)
zs_eff = z_recon * L_eff ./ (L_eff + z_recon);
d_window_sensor_to_sensor = zs_eff(1);
fprintf('Distances:\n   Sensor window to sensor: \t%.2f mm\n',d_window_sensor_to_sensor/mm);
if n_windows_estimation == 2
   d_window_source_to_sensor = zs_eff(3) + t2 - t2_eq;
   d_between_windows = d_window_source_to_sensor - t2 - d_window_sensor_to_sensor;
   d_window_source_to_source = L_eff - zs_eff(4);
   fprintf('   Window source to sensor: \t%.2f mm\n',d_window_source_to_sensor/mm);
   fprintf('   Between windows: \t\t%.2f mm\n',d_between_windows/mm);
   fprintf('   Window source to source: \t%.2f mm\n',d_window_source_to_source/mm);
end

% Magnification range
m_min = L_eff/(L_eff - zs_eff(2));
if n_windows_estimation == 2
   m_max = L_eff/(L_eff - zs_eff(3));
else % Guess that the window is right after the source
   zs_eff(3) = L_eff - t1eff; 
   m_max = L_eff/t1eff; 
end
fprintf('The magnification range in the zone between windows goes from %.1f to %.1f (source to sensor)\n',m_max,m_min)

% Diagonal of area under test
diag_side_sensor = diag_sensor / m_min;
diag_side_source = diag_sensor / m_max;
fprintf('The diagonal of the imaged zone goes from %.1f mm to %.1f mm.\n',diag_side_source/mm,diag_side_sensor/mm)

% Resolution range
NAs = sin(atan(diag_sensor./zs_eff(2:3)/2));
resol_sensor = max(1.22*wavel/NAs(1), dx_sensor/m_min);
resol_source = max(1.22*wavel/NAs(2), dx_sensor/m_max);
syms ds z w L pix
z_best_resol = solve(1.22*w/sin(atan(ds/z/2))==pix*(L-z)/L,z);
z_best_resol = min(double(subs(z_best_resol,[w ds pix L],[wavel diag_sensor dx_sensor L_eff])));
NA_best = sin(atan(diag_sensor./z_best_resol/2));
best_resol = [1.22*wavel/NA_best, dx_sensor*(L_eff-z_best_resol)/L_eff];
fprintf('The resolution goes from %.1f um to %.1f um, having a minimum of %.2f um at z=%.2f mm.\n',resol_source/um,resol_sensor/um,best_resol(1)/um,z_best_resol/mm)
