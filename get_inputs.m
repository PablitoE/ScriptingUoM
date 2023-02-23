% This function reads the entries from a text file getting the input
% arguments for the processing of an hologram.

function argument_str = get_inputs(filename)
% Default values
argument_str.image_path = '';    % String
argument_str.bgnd_path = [];     % Empty or string
argument_str.output_path = './';
argument_str.zeval = [];         % Evaluable string or whitespace separated values
argument_str.downsampling = 1;
argument_str.crop_factor = 1;
argument_str.z_autofoc_ini = []; % Linear range of analysis
argument_str.z_autofoc_end = 0;
argument_str.z_autofoc_n = 1;
argument_str.z_span = 0;         % Range of analysis around zs of interest, only used if z_autofoc_ini is empty
argument_str.z_resol = 0;
argument_str.z_interest = 1;     % Evaluable string
argument_str.L = 1;
argument_str.dx = 1e-6;
argument_str.dy = 1e-6;
argument_str.wavelength = .355e-6;
argument_str.make_movie = false;
argument_str.n1 = 1;
argument_str.n2 = 1;
argument_str.t1 = 0;
argument_str.t2 = 0;
argument_str.size_windows = 40; % May be 2 inputs of the minimum and maximum sizes of windows in pixels
argument_str.sub_mean_im = true;
argument_str.auto_min_max_video = true;
argument_str.use_equiv_wavel = false;

fileID = fopen(filename);
while ~feof(fileID)
    tline = fgetl(fileID);
    tline = split(tline);
    switch tline{1}
       case 'image_path'
          if length(tline) > 1
            argument_str.image_path = join(string(tline(2:end)));
          end
       case 'bgnd_path'
          if length(tline) > 1
             argument_str.bgnd_path = join(string(tline(2:end)));            
             if strcmp(argument_str.bgnd_path,'[]')
                argument_str.bgnd_path = [];
             end
          end
       case 'output_path'
          if length(tline) > 1
             argument_str.output_path = join(string(tline(2:end)));
          end
       case 'zeval'
          if length(tline) > 2
            eval(['argument_str.zeval = [' strjoin(tline(2:end)) '];'])
          elseif length(tline) == 2
            eval(['argument_str.zeval = ' tline{2} ';'])
          else
             warning('Empty evaluation zs.')
          end
       case 'downsampling'
          argument_str.downsampling = str2double(tline{2});
       case 'crop_factor'
          argument_str.crop_factor = str2double(tline{2});
       case 'z_autofoc'
          if strcmp(tline{2},'[]')
             argument_str.z_autofoc_ini = [];
          else
             argument_str.z_autofoc_ini = str2double(tline{2});
          end
          argument_str.z_autofoc_end = str2double(tline{3});
          argument_str.z_autofoc_n = str2double(tline{4});
       case 'z_interest'
          argument_str.z_span = str2double(tline{2});
          argument_str.z_resol = str2double(tline{3});
          argument_str.z_interest = eval(tline{4});
       case 'L'
          argument_str.L = str2double(tline{2});
       case 'dx'
          argument_str.dx = str2double(tline{2});
       case 'dy'
          argument_str.dy = str2double(tline{2});
       case 'wavelength'
          argument_str.wavelength = str2double(tline{2});
       case 'make_movie'
          argument_str.make_movie = eval(tline{2});
       case 'windows'
          argument_str.n1 = str2double(tline{2});
          argument_str.t1 = str2double(tline{3});
          argument_str.n2 = str2double(tline{4});
          argument_str.t2 = str2double(tline{5});
       case 'size_windows'
          if length(tline) == 2
              argument_str.size_windows = str2double(tline{2});
          elseif length(tline) == 3
              argument_str.size_windows = str2double(tline(2:3));
          end
       case 'sub_mean_im'
          argument_str.sub_mean_im = eval(tline{2});
       case 'auto_min_max_video'
          argument_str.auto_min_max_video = eval(tline{2});
       case 'use_equiv_wavel'
          argument_str.use_equiv_wavel = eval(tline{2});
    end
end
fclose(fileID);