% Finds the frequency of the peak with subpixel resolution in FFT2. The
% result is between 0 and sqrt(2), where 1 corresponds to sampling_frequency/2.
% guard : minimum frequency
% nmaxs : amount of maximae to consider when finding peak of main frequency
% n_patch : number of pixels to consider besides peak when using LS to fit
%           a parabola
% do_plot : show FFT
% only_x : 1D processing
% do_Fienup : refine result using Fienup subpixel cross correlation

function [freq, err] = get_main_frequency(im, guard, nmaxs, n_patch, do_plot, only_x, do_Fienup)
if ~exist('guard','var') || isempty(guard)
   guard = 0.005;
end
if ~exist('nmaxs','var') || isempty(nmaxs)
   nmaxs = floor(.0001*numel(im));    % 0.01% of pixels
end
if ~exist('n_patch','var') || isempty(n_patch)
   n_patch = 1;
end
if ~exist('do_plot','var') || isempty(do_plot)
   do_plot = false;
end
if ~exist('only_x','var') || isempty(only_x)
   only_x = false;
end
if ~exist('do_Fienup','var')
   do_Fienup = false; % It is not working
end

debug = false;

[M,N] = size(im);
middle = [floor(M/2+1) floor(N/2+1)];
fim = fftshift(fft2(im));
% Filter out zero order
zo_r = floor(M*(1-guard)/2+1):ceil(M*(1+guard)/2+1);
zo_c = floor(N*(1-guard)/2+1):ceil(N*(1+guard)/2+1);
fim(zo_r,zo_c) = 0;

fim = abs(fim);
if only_x
   fim1d = mean(fim);
   fim1d(zo_c) = 0;
   [~,ic] = max(fim1d);
   ir = middle(1);
else
   % Find a number of maximae coarsly
   [~, inds] = sort(fim(:),'descend');
   % Get the vector that gets the minimum error to a grid of points made
   % from that vector
   [rs, cs] = ind2sub([M,N], inds(1:nmaxs));
   rs = rs - middle(1);
   cs = cs - middle(2);
   xs = cs / N;
   ys = rs / M;
   % Base grid
   Ng = fix(sqrt(nmaxs)/2)*2 + 2;
   [grid_x, grid_y] = meshgrid(1:Ng,1:Ng);
   grid_x = grid_x - ceil(Ng/2);
   grid_y = grid_y - ceil(Ng/2);
   grid_p = [grid_x(:) grid_y(:)]';
   errors = zeros(nmaxs, 2);
   in_grid = zeros(nmaxs,1, 'logical');
   for k_peak = 1:nmaxs
      % Generate grid using this peak
      peak = [xs(k_peak), ys(k_peak)];      
      % Rotate ans scale grid with matrix formed from peak
      theta = atan(peak(2) / peak(1));
      R = norm(peak) * [cos(theta) -sin(theta); sin(theta) cos(theta)];
      this_grid = R * grid_p;
      % Find distances to grid points of all extrema
      ds = sqrt((xs - this_grid(1, :)).^2 + (ys - this_grid(2, :)).^2);
      % If this is a good extrema, then, there are 4 points with much less
      % error than others. Then, I need to keep the good point with minimum
      % norm
      minds = min(ds, [], 2);
      [sorted_ds, ind_ds] = sort(minds);
      in_grid(k_peak) = sum(sorted_ds(1:4).^2) < sum(sorted_ds(5:8).^2)/10;
      errors(k_peak, 1) = sum(minds);
      errors(k_peak, 2) = sum(minds.^2);
      
      % Debugging plots, see the grids and points
      if debug
         figure,
         plot(xs, ys, 'bx'), hold on
         plot(peak(1), peak(2), 'go')
         plot(this_grid(1,:), this_grid(2,:), 'r.')
         title(sprintf('Error %f', errors(k_peak, 1)))
         a = 1;
      end
   end
   % Min error [L1, L2] 
   [~, extr_L1] = min(errors(:, 1));
   [err, extr_L2] = min(errors(:, 2));
%    if extr_L1 ~= extr_L2
%       warning('The estimation of period would be different when using different norms.')
%    end
   [ir, ic] = ind2sub([M,N], inds(extr_L1));
%    [ir, ic] = ind2sub([M,N], inds(extr_L2));
   % % Get the distance to the middle from a set of nmaxs extrema   
   % dist_to_middle = vecnorm([rs cs],2,2);
   % [~, ind] = min(dist_to_middle);
   % [ir, ic] = ind2sub([M,N], inds(ind));
end

if do_plot
   figure
   imagesc(log(fim)), hold on
   plot(ic,ir,'o')
   [rsd, csd] = ind2sub([M,N], inds(1:nmaxs));
   plot(csd(in_grid), rsd(in_grid), 'b.')
   plot(csd(~in_grid), rsd(~in_grid), 'r.')
end

% Get a patch to fit a parabole
patch = fim(ir-n_patch:ir+n_patch,ic-n_patch:ic+n_patch);
is = -n_patch:n_patch;
[ix,iy] = meshgrid(is,is);
ix = ix(:); iy = iy(:);
A = [ix.^2 iy.^2 ix.*iy ix iy ones(numel(ix),1)];
% Estimate in log values
logpatch = log(patch(:));
% LS
p = (A'*A)\(A'*logpatch);

freq_x = (2*p(2)*p(4) - p(5)*p(3))/(p(3)^2-4*p(1)*p(2)) + ic - floor(N/2 + 1);
freq_y = (2*p(1)*p(5) - p(4)*p(3))/(p(3)^2-4*p(1)*p(2)) + ir - floor(M/2 + 1);
if do_plot
   plot(freq_x + middle(2),freq_y + middle(1),'x')
end
freq_x = freq_x /(N/2); freq_y = freq_y /(M/2);

% Refine by subpixel cross correlation to shifted image
if do_Fienup
   % Fourier (not fftshifted) of image
   F_im = fft2(im);
   % Shifted version
   period = 1 / sqrt((freq_x/2)^2 + (freq_y/2)^2);
   theta = atan(freq_y / freq_x);
   deltar = period * sin(theta);
   deltac = period * cos(theta);
   Nr = ifftshift((-fix(M/2):ceil(M/2)-1));
   Nc = ifftshift((-fix(N/2):ceil(N/2)-1));
   [Nc,Nr] = meshgrid(Nc,Nr);
   F_shifted_im = F_im.*exp(1i*2*pi*(deltar*Nr/M+deltac*Nc/N));
   
   % Fienup subpixel registration within 0.01 pixels by specifying an upsampling parameter of 100
   usfac = 100;
   [output, Greg] = dftregistration(F_im, F_shifted_im, usfac);
   freq_dy = M/2 / output(3);
   freq_dx = N/2 / output(4);
   freq_y = freq_y + freq_dy;
   freq_x = freq_x + freq_dx;
end

freq = sqrt(freq_x^2 + freq_y^2);