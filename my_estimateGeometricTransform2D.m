function [tform, good_matches, err_struct] = my_estimateGeometricTransform2D(uv,xy,transformationType,max_err,verbose)
% xy, uv : mx2
% A'A = lambda^2 I
% tform is the transformation required to take uv back to xy
err_struct = struct;
if size(xy,1) < 3
   tform = NaN;
   good_matches = NaN;
   err_struct.median_all = NaN;
   err_struct.mean_good = NaN;
   err_struct.max_good = NaN;
   err_struct.n_good = NaN;
   err_struct.n_all = NaN;
   err_struct.std_good = NaN;
   return
end

if ~exist('verbose','var') || ~islogical(verbose)
   verbose = true;
end
do_max_err = exist('max_err','var') && isnumeric(max_err); % maximum tolerated error to select outliers
MAX_ITER = 10;
% Initial removal of bad matches based on distances from xy to uv
distances = sqrt(sum((xy - uv).^2,2));
good_matches = find(~isoutlier(distances));
if length(good_matches) < 3
   good_matches = (1:size(xy,1))';
end
delta_good_ones = true;
iter = 1;
try
while delta_good_ones && iter <= MAX_ITER
   tform = fitgeotrans(uv(good_matches,:),xy(good_matches,:),transformationType);
   [expected_u, expected_v] = transformPointsInverse(tform,xy(:,1),xy(:,2));
   expected_uv = [expected_u, expected_v];
   errors = sqrt(sum((expected_uv - uv).^2,2));
   if do_max_err
      new_good_matches = find(errors < max_err);
   else
      new_good_matches = find(~isoutlier(errors));
   end
   if length(new_good_matches) < 3
      [~,ind_err] = sort(errors);
      new_good_matches = sort(ind_err(1:3)); % Keep best points
   end
   delta_good_ones = (numel(new_good_matches) ~= numel(good_matches)) ||...
      any(new_good_matches ~= good_matches);
   good_matches = new_good_matches;
   iter = iter + 1;
end
if iter > MAX_ITER && delta_good_ones && verbose
   warning('Estimation of Geometric transform ended as reached the maximum number of allowed iterations.')
end
catch
   tform = NaN;
   good_matches = NaN;
   err_struct.median_all = NaN;
   err_struct.mean_good = NaN;
   err_struct.max_good = NaN;
   err_struct.n_good = NaN;
   err_struct.n_all = NaN;
   err_struct.std_good = NaN;
   return
end
err_struct.median_all = median(errors);
err_struct.mean_good = mean(errors(good_matches));
err_struct.max_good = max(errors(good_matches));
err_struct.n_good = numel(good_matches);
err_struct.n_all = numel(errors);
err_struct.std_good = std(errors(good_matches));
if verbose
   fprintf('Median error of all matches: %.2f\n', err_struct.median_all)
   fprintf('Mean error of good matches: %.2f\n', err_struct.mean_good)
   fprintf('Maximum error of good matches: %.2f\n', err_struct.max_good)
end