function [gr_ini, gr_end, im1_good] = detect_good_row_start_end(im1,im2,fill_rows)

WIDTH_MV = 5; % width of window for getting the mean value of fringes or neighbors

if nargin == 1
   im2 = [];
end
if nargin < 3
   fill_rows = false;   % By default, not filling bad rows
end
[M,N] = size(im1);
sep = floor(M/2);
im1_rv = mean(im1,2);
m1 = mean(im1_rv);
vstd = std(im1_rv);
im1_r = im1_rv < m1-3*vstd; % Detect very low intensity
im1_r_1 = im1_r(1:sep);
im1_r_2 = im1_r(sep+1:end);
% Initial good row
beginnings = strfind([0 0 im1_r'],[0 0 1 1]);
endings = strfind([im1_r' 0 0],[1 1 0 0]) + 1;
Nbad_fringes = length(beginnings);
if Nbad_fringes ~= length(endings)
   error('Bad detection of beggining and ending of bad rows.')
end
% Detect if bad fringes are in the upper and lower halves of the image
is_upper_bad_fringe = beginnings < sep;
Nbad_upper = sum(is_upper_bad_fringe);
Nbad_lower = Nbad_fringes - Nbad_upper;
% Use inpainting if fill_rows == true, and bad rows were detected
if fill_rows && ~isempty(beginnings) 
%    im1_good = im1;
%    mean_value_before = mean(im1_rv(max(beginnings-WIDTH_MV,1):beginnings-1));
%    mean_value_after = mean(im1_rv(endings+1:min(endings+WIDTH_MV,length(im1_r_1))));
%    im1_good(beginnings:endings,:) = (mean_value_before + mean_value_after)/2;
   mask = zeros(size(im1),'logical');
   mask(beginnings:endings,:) = true;
   im1_good = inpaintCoherent(im1,mask);
   ini_gr1 = [];
   
elseif Nbad_upper > 1 && Nbad_upper == Nbad_fringes % Is it possible to correct upper bad fringes? (No lower bad fringes)
   im1_good = im1;
   % Each fringe and space between fringes should be multiplied by a factor
   % that would reduce the error between the mean value of itself and neighbors
   mean_values = zeros(1,2*Nbad_fringes+1);
   mean_values(1) = mean(im1_rv(max(beginnings(1)-WIDTH_MV,1):beginnings(1)-1));
   mean_values(end) = mean(im1_rv(endings(end)+1:min(endings(end)+WIDTH_MV,length(im1_r_1))));
   for kfringe = 1:Nbad_fringes
      % Get mean value in fringes
      mean_values(2*kfringe) = mean(im1_rv(beginnings(kfringe):endings(kfringe)));
      % Get mean value between fringes
      if kfringe < Nbad_fringes
         mean_values(2*kfringe+1) = mean(im1_rv(endings(kfringe)+1:beginnings(kfringe+1)-1));
      end
   end
   % Cost function : (Delta(mv .* [1 a 1]))^2
   % Linear equations system : Gradient = 0
   A = 2*eye(2*Nbad_fringes-1) - diag(ones(2*Nbad_fringes-2,1),-1) - diag(ones(2*Nbad_fringes-2,1),1);
   A = A .* mean_values(2:end-1);
   b = zeros(2*Nbad_fringes-1,1);
   b(1) = mean_values(1);
   b(end) = mean_values(end);
   a = A\b;
   for kfringe = 1:Nbad_fringes
      % Correct fringe
      im1_good(beginnings(kfringe):endings(kfringe),:) = im1(beginnings(kfringe):endings(kfringe),:) * a(kfringe*2-1);
      % Correct space between fringes
      if kfringe < Nbad_fringes
         im1_good(endings(kfringe)+1:beginnings(kfringe+1)-1,:) = im1(endings(kfringe)+1:beginnings(kfringe+1)-1,:) * a(kfringe*2);
      end
   end
   ini_gr1 = [];
else
   im1_good = [];
   ini_gr1 = find(im1_r_1,1,'last'); % Find bad rows
end
if isempty(ini_gr1)
   gr1_ini = 1;
else
   im1_r_1(1:ini_gr1) = true; % fill good portion before bad rows
   gr1_ini = find(~im1_r_1);   % find first good row
end
% Ending good row
gr1_end = find(im1_r_2,1)-1 + sep;   % find first bad row
if isempty(gr1_end) || fill_rows
   gr1_end = M;
end

if ~isempty(im2)
   im2_r = mean(im2,2);
   m2 = mean(im2_r);
   im2_r = im2_r < m2/3;
   im2_r_1 = im2_r(1:sep,:);
   im2_r_2 = im2_r(sep+1:end);
   % Initial good row
   ini_gr2 = find(im2_r_1,1);
   if isempty(ini_gr2)
      gr2_ini = 1;
   else
      im2_r_1(1:ini_gr2) = true; % fill good portion before bad rows
      gr2_ini = find(~im2_r_1);   % find first good row
   end
   % Ending good row
   gr2_end = find(im2_r_2,1)-1 + sep;   % find first bad row
   if isempty(gr2_end)
      gr2_end = M;
   end
   gr_ini = max(gr1_ini,gr2_ini);
   gr_end = min(gr1_end,gr2_end);
else
   gr_ini = gr1_ini;
   gr_end = gr1_end;
end