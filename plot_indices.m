% Function for analising the results in mat file obtained by Single_ASM_j
% in combination with join_Indices when using job array
%
% filename : path and filename of MAT file with found indices
% zs : z reconstruction positions for the indices
% size_wins : size of windows for getting each column of indices
% select_size : possible to choose some of the windows to plot, 'sum', or
%               all [].
% lp : [] or (0,1) for lowpass filter

function plot_indices(filename,zs,size_wins,select_size,lp, ind_names)
close all
if ~exist(filename,'file')
   error('The input file name was not found.')
end
if ~strcmp(filename(end-3:end),'.mat')
   error('The input file is not a MAT file.')
end
variableInfo = who('-file',filename);
if ismember('all_indices',variableInfo)
   load(filename,'all_indices');
elseif ismember('joined_indices',variableInfo)
   load(filename,'joined_indices');
   all_indices = joined_indices;
   clear joined_indices
else
   error('Wrong file.')
end
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
elseif strcmp(select_size,'sum')
   size_wins = ['Sum of ' num2str(size_wins)];
end

fig_indices = figure(1); 
fig_indices_opos_maxmin = figure(2);
if nargin < 6 || isempty(ind_names)
   ind_names = {'Entropy', 'L1', 'Gradient', 'Laplacian', 'Variance', 'Tamura', 'CC', 'Gini'}; % , 'RC'
end

for ki = 1:length(all_indices)
   mat_ind = all_indices{ki};
   [Nz,Nw,Nq] = size(mat_ind);
   loc_nan = isnan(mat_ind) | isinf(mat_ind);
   mat_ind(loc_nan) = 0;
   if isnumeric(select_size)
      mat_ind = mat_ind(:,select_size,:);
   end
   if ~isempty(lp)
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
      for kq = 1:Nq
         try
            mat_ind(:,:,kq) = lowpass(mat_ind(:,:,kq),lp);
         catch
            mat_ind(:,:,kq) = filtfilt(lpFilt,mat_ind(:,:,kq));
         end
      end
   end
   if strcmp(select_size,'sum')
      mat_ind = sum(mat_ind,2);
   end
   figure(fig_indices)   
   ki_sp = mod(ki-1,8)+1;
   subplot(2,4,ki_sp)
   hold on
   plot(zs, mat_ind(:,:,1))
   title([get(get(gca,'title'),'String') ' ' ind_names{ki}])
   revise_last_legend(size_wins) % Get proper legend
   figure(fig_indices_opos_maxmin)
   subplot(2,4,ki_sp)
   hold on
   plot(zs, mat_ind(:,:,2))
   title([get(get(gca,'title'),'String') ' ' ind_names{ki}])
   revise_last_legend(size_wins) % Get proper legend
end
% saveas(fig_indices,'.\results\Indices.fig')
end

function revise_last_legend(v)
   leg_obj = get(gca,'legend');
   leg_htal_array = get(leg_obj,'String');
   good_leg = num2str(v(:));
   if isempty(leg_htal_array)
      legend(good_leg)
   else
      if isnumeric(v)
         N = length(v);
         for k = 1:N
            leg_htal_array{end-N+k} = good_leg(k,:);
         end      
      elseif ischar(v)
         leg_htal_array{end} = good_leg;
      end         
      legend(leg_htal_array)
   end
end