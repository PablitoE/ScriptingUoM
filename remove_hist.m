% Remove histogram that the software of Photron may add to images
% path_ims = 'D:\HoloICE\Holograms\Pollen\Photron\ash_inv3Ms_maxpower_01\Calibration';
path_ims = 'D:\HoloICE\Holograms\Navitar_test_Photron\No tubes\USAF 1951';
format = 'png';

list_dir = dir(path_ims);
n_ims = 0;
names_ims = {};
for k=1:length(list_dir)
   file = list_dir(k);
   try
      fullpath_file = fullfile(path_ims,file.name);
      warning('off','all')
      a = imfinfo(fullpath_file);
      warning('on','all')
      n_ims = n_ims + 1;
      names_ims{n_ims} = file.name;
      im = imread(fullpath_file);
      
      non_grey_color = range(im, 3) > 0;
      fonts_and_lines = all(im == 255 | im == 128, 3);
      mask_inpaint = non_grey_color | fonts_and_lines;
      clear non_grey_color fonts_and_lines
      % Histogram only in lower left corner
      [M, N] = size(mask_inpaint);
      mask_inpaint(1:M-fix(M/3), :) = 0;
      mask_inpaint(:, fix(N/3):N) = 0;
      im = im(:, : , 1);
      im_inpaint = inpaintCoherent(im, mask_inpaint);
      
      % Save clean image
      clean_file = fullfile(path_ims, sprintf('%s_noHist.png', file.name(1:end-4)));
      imwrite(im_inpaint, clean_file)
   catch e
        fprintf('The identifier was:\n%s',e.identifier);
        fprintf('There was an error! The message was:\n%s',e.message);
   end
end