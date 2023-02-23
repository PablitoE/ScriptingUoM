% The input FT_Holo must have been fftshifted

function fieldOut = easy_prop(k, zeval, root, FT_Holo, padding,verbose)
if ~exist('verbose','var')
   verbose = true;
end

[Nr, Nc] = size(FT_Holo);
if verbose, fprintf('Prop ASM : '), end

% Propagation to desired z
phaseFilter = exp(1j*k * zeval * root);
if verbose, fprintf('phase filter ready... '), end
FT_Holo = FT_Holo .* phaseFilter;
clear phaseFilter
if verbose, fprintf('mult ready... '), end

% Take the ifft2 which completes the transform
fieldOut = ifft2(ifftshift(FT_Holo));
if verbose, fprintf('field ready... '), end
fieldOut = single(fieldOut); % it is single already

% Crop due to padding
fieldOut = fieldOut(padding(1)+1:Nr-padding(1), padding(2)+1:Nc-padding(2));
if verbose, fprintf('Done.\n'), end
