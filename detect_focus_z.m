% Function for analising the results in mat file obtained by
% Focus_detection_single_plane
%
% indices : obtained focus detection indices
% zs : z reconstruction positions for the indices
% size_wins : size of windows for getting each column of indices
% select_size : possible to choose some of the windows or all [] to sum results.
% lp : [] or (0,1) for lowpass filter

function z_focus = detect_focus_z(all_indices,zs,size_wins,select_size,lp,do_plot,out_folder,name_image)
close all
if nargin < 4 || isempty(select_size)
   select_size = 1:length(size_wins);
end
if nargin < 5 || ~isnumeric(lp) || isempty(lp)
   lp = [];
else
   lpFilt = designfilt('lowpassiir','FilterOrder',8, ...
         'PassbandFrequency',lp,'PassbandRipple',0.2, ...
         'SampleRate',1);
end
if isnumeric(select_size)
   size_wins = size_wins(select_size);
end
if nargin < 6 || isempty(do_plot)
   do_plot = false;
end
if nargin < 7 || isempty(out_folder)
   out_folder = '.\';
end
if nargin < 8 || isempty(name_image)
   name_image = '';
end
if do_plot
   fh = figure('units','normalized','outerposition',[0 0 1 1]);
end
Nind = length(all_indices);
the_maxs = zeros(Nind, 1);
z_position = zeros(Nind, 1);
Fs = cell(Nind, 1);
for ki = 1:Nind
   mat_ind = all_indices{ki};
   if isnumeric(select_size)
      mat_ind = mat_ind(:,select_size,:);
   end   
   [~,Nw,Nq] = size(mat_ind);
   loc_nan = isnan(mat_ind) | isinf(mat_ind);
   mat_ind(loc_nan) = 0;
   % If non uniform sampling, resample
   dzs = diff(zs);
   if any(dzs ~= dzs(1))
      mindz = min(dzs(dzs > 0));
      new_zs = zs(1):mindz:zs(end);
      new_mat_ind = zeros(length(new_zs),Nw,Nq);
      for kq = 1:Nq
         new_mat_ind(:,:,kq) = interp1(zs,mat_ind(:,:,kq),new_zs);
      end
   end
   if ~isempty(lp)      
      for kq = 1:Nq
         try
            new_mat_ind(:,:,kq) = lowpass(new_mat_ind(:,:,kq),lp);
         catch
            new_mat_ind(:,:,kq) = filtfilt(lpFilt,new_mat_ind(:,:,kq));
         end
      end
   end
   
   mat_ind = sum(new_mat_ind,2);    % Sum results from different window sizes
   % Remove trend (usually descending)
   mat_ind_det = detrend(mat_ind(:, 2)); % Don't use maximums, use sum(x^EXP_MAX_MIN)   
   [~, kz] = max(mat_ind_det);
   % Refine the maximum position
   F = griddedInterpolant(new_zs,mat_ind_det,'spline');
   Fs{ki} = F;
   if kz == 1
      z_position(ki) = new_zs(1);
   elseif kz == length(mat_ind_det)
      z_position(ki) = new_zs(end);
   else
      z_position(ki) = fminbnd(@(x)-F(x),new_zs(max(kz-2,1)), new_zs(min(kz+2,numel(new_zs))));
%       z_position(ki) = fminsearch(@(x)-F(x),new_zs(kz));
%       if z_position(ki) < new_zs(1) || z_position(ki) > new_zs(end)
%          z_position(ki) = fminbnd(@(x)-F(x),new_zs(1), new_zs(end));
%       end
   end
   the_maxs(ki) = F(z_position(ki));

   if do_plot
      subplot(1,Nind,ki)
      plot(new_zs, mat_ind(:, 2),':'), hold on     % No detrending
      plot(new_zs, mat_ind_det)                    % Detrended
      plot(z_position(ki),the_maxs(ki),'o')        % Maximum of this index
      plot([1 1]*z_position(ki),ylim,'--')         % Add vertical line there
   end
end
% Weighted mean
weigths = the_maxs / sum(the_maxs);
z_focus = sum(z_position.*weigths);
% Weigh the amplitude of other indices in every max
all_weights = zeros(Nind);
for ki = 1:Nind
   the_F = Fs{ki};
   all_weights(ki, :) = the_F(z_position);
end
all_weights = all_weights - min(all_weights(:));
all_weights = all_weights / sum(sum(all_weights));
z_focus2 = sum(sum(all_weights .* z_position'));
if do_plot
   for ki = 1:Nind
      subplot(1,Nind,ki)
      plot([1 1]*z_focus,ylim,':g')       % Vertical line at definitive focus point
      plot([1 1]*z_focus2,ylim,':m')       % Vertical line at 2nd possible focus point
   end
   name_file_fig = sprintf('focus_detection_%s_%5.0fmm.jpg',name_image,z_focus*1e3);
   saveas(fh,fullfile(out_folder,name_file_fig))
   fprintf('Figure file created: %s', name_file_fig)
   close(fh)
end
end