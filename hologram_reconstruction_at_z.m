function Psi = hologram_reconstruction_at_z(hologram, z, wavelength, dx)
[Nr, Nc] = size(hologram);
hologram = double(hologram) - mean(hologram(:));
padding = [Nr/2, Nc/2];
hologram = padarray(hologram, padding);
Nr = Nr * 2; Nc = Nc * 2;
verbose = false;
k = 2*pi / wavelength;

FT_Holo = fftshift(fft2(hologram));
clear hologram
root = get_ASM_root(dx,dx,Nc,Nr,wavelength);
Psi = easy_prop(k, z, root, FT_Holo, padding, verbose);
end