% Function for getting the period of a periodic image (X and Y) by getting
% the maximum peak after the zero order peak.
% Input 'guard' is a number from 0 to 1 representing the high-pass filter
% cutoff.
function period = im_peak_period(im, dx, guard)
[M, N] = size(im);
if M~=N
   error('The image must be square.')
end
if ~exist('guard','var')
   guard = 0.05;
end

FT = fft2(im);
FT = FT(:, 1:floor(N/2));
mFT = mean(abs(FT));
M = length(mFT);
mFT(1:min(ceil(guard*M),M)) = 0;  % Keep positive frequencies
[~, ipeak] = max(mFT);
df = 1/dx/N;
fpeak = (ipeak-1)*df;
period = 1/fpeak;