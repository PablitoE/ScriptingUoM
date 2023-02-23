% Get the z values of a divergent beam that would obtain a constant
% step-to-resolution ratio.
% start : initial z
% stop : last z
% N : number of steps. Total of N+1 points from start to stop
% L : distance in m from sensor to source (equivalent in case of having
% windows in the middle)
% D : maximum dimension of the sensor
% wavelength
% Output: 
%  zs : locations in z with constant K
%  K : step-to-resolution ratio

function [zs, K] = sample_z_sph(start, stop, N, L, D, wavelength)
% Solve the problem in planar wavefront and then go back to divergent
start = z_divergent2planar(start, L);
stop = z_divergent2planar(stop, L);
% Solve K that gets f(K) == 0
fun = @(x) f_to_make_null(x, start, stop, N, D, wavelength);
% The seed is obtained by using logspace
% z_log = logspace(log10(start), log10(stop), N+1);
% first_step = z_log(2) - z_log(1);
% Kseed = first_step / get_resolution(start, L, D, wavelength);
Kseed = 1;
% options = optimoptions('fsolve', 'TolFun', .1);
options = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt', 'TolFun', 1e-9, 'TolX', 1e-9);
K = fsolve(fun, Kseed, options);
zs = get_zs_from_K(K, start, N, D, wavelength);
zs(end) = stop;
zs = z_planar2divergent(zs, L);
end

function out = f_to_make_null(K, start, stop, N, D, wavelength)
zs = get_zs_from_K(K, start, N, D, wavelength);
out = stop - zs(end);
end

function zs = get_zs_from_K(K, start, N, D, wavelength)
zs = zeros(N,1);
zs(1) = start;
for k=1:N-1
   delta = get_resolution(zs(k), D, wavelength);
   zs(k+1) = zs(k) + K*delta;
end
end

function dz = get_resolution(z, D, wavelength)
% Get resolution in a propagation of planar wave
NA = sin(atan(D/2/z));
dz = wavelength / NA^2;
end

function z_planar = z_divergent2planar(z_divergent,L)
   z_planar = z_divergent*L./(L-z_divergent);
end

function z_divergent = z_planar2divergent(z_planar,L)
   z_divergent = z_planar*L./(L+z_planar);
end