function [im_d, dx, dy, Nc, Nr] = downsample_image(im,downsampling,dx0,dy0)
[Nr0, Nc0] = size(im);
if exist('dx0','var') && ~isempty(dx0)
   dx = dx0 * downsampling;
else
   dx = [];
end
if exist('dy0','var') && ~isempty(dy0)
   dy = dy0 * downsampling;
else
   dy = [];
end
% im = im(1:downsampling:end, 1:downsampling:end); % Aliasing
% im = imresample([1, 1], im, [1,1]/downsampling,'cubic');  % Too much memory
F = griddedInterpolant(single(im));
cc = (1:downsampling:Nc0)';
rr = (1:downsampling:Nr0)';
im_d = F({rr,cc});
clear F
[Nr, Nc] = size(im_d); % [Ny, Nx]
% Resize the image to make it FFT friendly
if mod(Nc,2), im_d = im_d(:,1:end-1); Nc = Nc -1; end
if mod(Nr,2), im_d = im_d(1:end-1,:); Nr = Nr -1; end
