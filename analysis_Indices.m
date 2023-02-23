% Analysis of Indices.mat, results from Single_ASM_j. join_Indices was used
% to join the matrices coming from different nodes in JASMIN.
path_to_file = 'D:\HoloICE\Scripting\results\PhotronUSAFCali\Indices.mat';

load(path_to_file)

if ~isvarname('autof_Ws')
   mm = 1e-3;
   um = 1e-6;
   nm = 1e-9;

   z_autofoc = [10e-3 72e-3 4000]; % or []
   z_resol = []; z_span = []; z_interest = [];
   L = 80*mm;
   size_windows = [7, 160];

   dx0 = 2.2*um;
   dy0 = 2.2*um;
   Nr0 = 9176; 
   Nc0 = 13264;
   wavelength = 355*nm;

   n1 = 1.4761;
   n2 = 1.4761;
   t1 = 5*mm;
   t2 = 12*mm;
   % Dealing with windows effects, equivalent distance source-sensor
   t1_eq = t1 / n1;
   t2_eq = t2 / n2;
   L = L - t1 + t1_eq - t2 + t2_eq;
   wavelength_equivalent = wavelength;

   D = max(Nr0 * dy0, Nc0 * dx0);

   % Propagation distances
   if ~isempty(z_autofoc)
   %     z_autofoc = linspace(str_input.z_autofoc_ini,str_input.z_autofoc_end,str_input.z_autofoc_n);
       [z_autofoc, step_to_resolution_ratio] = sample_z_sph(z_autofoc(1),z_autofoc(2),z_autofoc(3),L,D,wavelength_equivalent);
       fprintf('The step-to-resolution ratio is %.2f.\n', step_to_resolution_ratio);
   elseif ~isempty(z_resol)
       z_autofoc = (-z_span/2:z_resol:z_span/2)' + z_interest;
       z_autofoc = z_autofoc(:);
   end
   if any(z_autofoc >= L)
      z_autofoc = z_autofoc(z_autofoc < L);
   end
end

% 1 Get a look of the results from using different sizes of windows (extrema and mean N-power)
% plot_indices(path_to_file,z_autofoc,autof_Ws)
% 2
% plot_indices(path_to_file,z_autofoc,autof_Ws,'sum')   % Sums results over all windows
% 3
plot_indices(path_to_file,z_autofoc,autof_Ws,'sum',0.2)   % Filtered