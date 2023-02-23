% Axial sampling with constant step to resolution ratio: K.
% A reconstruction based on a planar wavefront is considered.
% Taken from config.axial_sampling
function zs = axial_sampling(zMin, zMax, im_shape, pixel_size, K, wavelength)
   maxDimSensor = max(pixel_size.*im_shape);
   
   zs = z_step_to_resol(zMin,zMax,[],maxDimSensor,wavelength,[],K);
end

function [zs, K] = z_step_to_resol(start, stop, dz0, D, wavelength, N, K)
% Get the z values of a planar wave that would obtain a constant
% step-to-resolution ratio. All inputs in metres.
% start : initial z
% stop : last z
% dz0 : initial step
% D : maximum dimension of the sensor
% wavelength : laser wavelength
% [N] : number of z positions in output (dz0 is ignored)
% [K] : step-to-resolution ratio (dz0 and N are ignored)
%
% Output:
%  zs : locations in z with constant K
%  K : step-to-resolution ratio

case_sampling = 'dz0';
if exist('N','var') && ~isempty(N) && uint32(N)==N, case_sampling = 'N'; end
if exist('K','var') && ~isempty(K) && K>0, case_sampling = 'K'; end

switch case_sampling
   case 'dz0'
      K = dz0 / get_axial_resolution(start, D, wavelength);
      zs = get_zs_from_K_until_stop(K, start, D, wavelength, stop);
   case 'N'
      % Solve K that gets f(K) == 0
      fun = @(x) f_to_make_null(x, start, stop, N, D, wavelength);
      Kseed = 1;
      % options = optimoptions('fsolve', 'TolFun', .1);
      options = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt', 'TolFun', 1e-9, 'TolX', 1e-9,'Display','off');
      [K,fval,exitflag,output] = fsolve(fun, Kseed, options);
      if exitflag < 1
         error('The sampling solution is not reliable. Choose different sampling mode or inputs.')
      end
      zs = get_N_zs_from_K(K, start, N, D, wavelength);
      zs(end) = stop;
   case 'K'
      zs = get_zs_from_K_until_stop(K, start, D, wavelength, stop);
end
end

function zs = get_zs_from_K_until_stop(K, start, D, wavelength, stop)
zs = start;
while stop > zs(end)
   delta = get_axial_resolution(zs(end),D, wavelength);
   zs(end+1) = zs(end) + K*delta;
end
if length(zs) > 2
   zs = zs(1:end-1);
end
end

function out = f_to_make_null(K, start, stop, N, D, wavelength)
zs = get_N_zs_from_K(K, start, N, D, wavelength);
out = stop - zs(end);
end

function zs = get_N_zs_from_K(K, start, N, D, wavelength)
zs = zeros(N,1);
zs(1) = start;
for k=1:N-1
   delta = get_axial_resolution(zs(k), D, wavelength);
   zs(k+1) = zs(k) + K*delta;
end
end

function dz = get_axial_resolution(z, D, wavelength)
% Get resolution in a propagation of planar wave
NA = sin(atan(D/2/z));
dz = wavelength / NA^2;
end