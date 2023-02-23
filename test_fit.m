% try fit airy
the_x = [6, 6, 2, 0, 2, 2];
[cc, rr] = meshgrid(1:12, 1:13);
y = airymodel(the_x,[cc(:),rr(:)]);
im = reshape(y,[13,12]);
im_nans = im; im_nans(randlocs) = nan;
[mma,mna] = fit_im(im_nans,@airymodel);