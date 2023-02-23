% Script for making a movie of a particular portion of the FOV of a
% sequence of holograms
configfile = 'D:\HoloICE\Holograms\Rime splintering\20230216\dropsOnFrost_9X_04_C001H001S0002\local_win.cfg';
outputImageDir = 'D:\HoloICE\Holograms\Rime splintering\20230216\dropsOnFrost_9X_04_C001H001S0002\recon';
imageNamesFormat = 'dropsOnFrost_9X_04_C001H001S000200%04d.png';
cropCoords = [151, 600, 251, 550];  % [row_start, row_end, col_start, col_end]
reconstructedToZ = 9.55e-3;

outputImageFactorSize = 2;
sequenceNumber = 1;
outputImageSuffix = '_temp';
shouldOverwrite = false;
pixsize = 20e-6 / 9;
scaleValues = 1.3; % Increases amplitude multiplying by scaleValues

width_mm = ((cropCoords(4) - cropCoords(3) + 1) * pixsize) * 1e3;
height_mm = ((cropCoords(2) - cropCoords(1) + 1) * pixsize ) * 1e3;
fprintf("Size of images is : %f mm x %f mm\n", width_mm, height_mm)

pathsToTempImages = makeMovieFromHologramSequence(configfile, sequenceNumber, outputImageDir,...
   cropCoords, outputImageFactorSize, outputImageSuffix, reconstructedToZ, scaleValues, shouldOverwrite);

imageNamesFormat = strrep(imageNamesFormat,'.png', [outputImageSuffix '.png']);
commandFFMPEG = sprintf('ffmpeg -framerate 10 -i "%1$s%2$c%3$s" -c:v libx264 -pix_fmt yuv420p "%1$s%2$cout.mp4"', outputImageDir, filesep, imageNamesFormat);
out = system(commandFFMPEG);

if out == 0
   for ktemp = 1:length(pathsToTempImages)
      delete(pathsToTempImages{ktemp})
   end
end