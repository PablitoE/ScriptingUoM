% Join the outputs of job array result

function vout = join_Indices(folder, work_anyway)
if nargin == 1
   work_anyway = false; % Default: the amount of files must be equal to the maximum number in the names of the mat files.
end
fprintf('Looking for Indices_xxxx.mat files in %s.\n', folder);
listing = dir(folder);
if ~isempty(listing)
   Nf = length(listing);
   valid_files = 0;
   max_num = 0;
   found_nums = zeros(Nf,1,'logical');
   for kf = 1:Nf
      % If name satisfies structure of mat files...
      if ~listing(kf).isdir && strcmp(listing(kf).name(1:8),'Indices_') && strcmp(listing(kf).name(13:end),'.mat')
         num = str2double(listing(kf).name(9:12));
         found_nums(num) = true;    % Remember files that were found
         max_num = max(num,max_num);
         valid_files = valid_files + 1;
      end
   end
   found_nums = found_nums(1:max_num); % Remove empty flags
   % Case of not having the same amount of files as said by numbers in the
   % names of the file, and not working anyway. Raise error.
   if valid_files ~= max_num && ~work_anyway
      not_found = find(~found_nums);
      fprintf('Files not found : ')
      fprintf('%d ', not_found)
      fprintf('\n')
      save(fullfile(folder, 'not_found.mat'), 'not_found');
      error(['The numbers in the files names are not correct. %d valid files were found, whereas %d is the maximum number in the files. A file not_found.mat'...
         ' was created to be used by script run_jobs_errored using %d cores.'],valid_files,max_num,length(not_found))
   end
   % Get nums of valid files
   found_nums = find(found_nums);
   if numel(found_nums) ~= valid_files
      error('Something is wrong with the amount of valid files.')
   end
   % Load first valid file
   filename = fullfile(folder,sprintf('Indices_%04d.mat', found_nums(1)));
   loaded = load(filename);      
   joined_indices = loaded.all_indices;
   if valid_files > 1
      for kf = 2:valid_files
         filename = fullfile(folder,sprintf('Indices_%04d.mat',found_nums(kf)));
         loaded = load(filename);
         joined_indices = cellfun(@(x,y)[x;y],joined_indices,loaded.all_indices,'UniformOutput',false);      
      end
   end
%    save(fullfile(folder,'Indices.mat'),'joined_indices')
   % Save other information contained in all files.
   saving = loaded;
   saving = rmfield(saving,'all_indices');
   saving.joined_indices = joined_indices;
   save(fullfile(folder,'Indices.mat'),'-struct','saving')
   for kf = 1:valid_files
      filename = fullfile(folder,sprintf('Indices_%04d.mat',found_nums(kf)));
      delete(filename)
   end
else
   error('Empty directory.')
end
vout = 1;