%    Used to make a move of a hologram series or a series of .png files

%    Copyright (C) 2014 Jacob P. Fugal
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%    When using this program as part of research that results in
%    publication, acknowledgement is appreciated. The following citation is
%    appropriate to include:
%
%    Fugal, J. P., T. J. Schulz, and R. A. Shaw, 2009: Practical methods
%    for automated reconstruction and characterization of particles in
%    digital inline holograms, Meas. Sci. Technol., 20, 075501,
%    doi:10.1088/0957-0233/20/7/075501.
%
%    Funding for development of holoViewer at Michigan Tech (Houghton,
%    Michigan, USA) provided by the US National Science Foundation (US-NSF)
%    Graduate Research Fellowship Program, and NASA's Earth Science
%    Fellowship Program. Funding for development at the National Center for
%    Atmospheric Research (NCAR, Boulder, Colorado, USA) provided by the
%    US-NSF. Funding for development at the Max Planck Institute for
%    Chemistry (MPI-Chemistry, Mainz, Germany) and the Johannes Gutenberg
%    University of Mainz (Uni-Mainz, Mainz, Germany), provided by
%    MPI-Chemistry, the US-NSF, and the Deutsche Forschungsgesellschaft
%    (DFG, the German Science Foundation).
%
%    Please address questions or bug reports to Jacob Fugal or to Matthew
%    Beals at fugalscientific (at) gmail (dot) com or mjbeals (at) mtu
%    (dot) edu respectively

function allnewpaths = makeMovieFromHologramSequence(configfile, sequenceNumber, outputImageDir,...
   cropCoords, outputImageFactorSize, outputImageSuffix, reconstructedToZ, scaleMax, shouldOverwrite)

% Set defaults if not specified.
if ischar(configfile), cfg = config(configfile); end
if ~exist('sequenceNumber','var') || isempty(sequenceNumber), sequenceNumber = 1; end
if ~exist('outputImageDir','var') || isempty(outputImageDir), outputImageDir = '.'; end
if ~exist('cropCoords','var') || isempty(cropCoords), cropCoords = []; end
if ~exist('outputImageFactorSize','var') || isempty(outputImageFactorSize), outputImageFactorSize = 1; end
if ~exist('outputImageSuffix','var') || isempty(outputImageSuffix), outputImageSuffix = '_'; end
if ~exist('reconstructedToZ','var') || isempty(reconstructedToZ), reconstructedToZ = 0; end
if ~exist('scaleMax','var') || isempty(scaleMax), scaleMax = 1; end
if ~exist('shouldOverwrite','var') || isempty(shouldOverwrite), shouldOverwrite = false; end

% Check for the existence of the output directory
if ~exist(outputImageDir, 'dir'), mkdir(outputImageDir); end

% Get the hologram list
holoFilenames = cfg.seq_list(sequenceNumber);

% Load the first hologram to get the image size
hologram = imread([cfg.path filesep holoFilenames{1}]);
[Ny, Nx] = size(hologram);

% Begin progress monitor
%uS = etd(clock, 1, length(holoFilenames), 60);
configfile = cfg.config_path;
allnewpaths = cell(length(holoFilenames), 1);

parfor cnt=1:length(holoFilenames)
   % Make the save filename and check for existence
   stime = clock;
   saveFilename = config.replaceEnding(holoFilenames{cnt}, '.png', [outputImageSuffix '.png']);
   allnewpaths{cnt} = [outputImageDir filesep saveFilename];
   if ~shouldOverwrite 
      if exist(allnewpaths{cnt},'file')
         continue;
      end
   end
   
   % Get the hologram we will display
   %    cfg.current_holo = holoFilenames{cnt};
   % Do the prefilters
   cfg = config(configfile);
   cfg.current_holo = holoFilenames{cnt};
   preimage = img;
   preimage.filter_handle = 'prefilters';
   preimage.config_handle = cfg;
   preimage.raw_image = cfg.full_holo_path;
   hologram = preimage.ampFiltered;
   
   % Reconstruct to some plane if requested to:
   if reconstructedToZ ~= 0
      pEngine = Propagator;
      pEngine.config_handle = cfg;
      pEngine.preconstruct(hologram);
      pEngine.unregisterListeners;
      slice    = pEngine.slice(reconstructedToZ);
   else
      slice = hologram;
   end
   % Apply postfilters
   postimage = img;
   postimage.filter_handle = 'postfilters';
   postimage.config_handle = cfg;
   postimage.unregisterListeners;
   postimage.raw_image = slice;
   postim   = postimage.ampEnhanced;   
   postim = scaleMax * postim;
   % Convert from double from [0 1] to uint8 from [0 255]
   saveim   = im2uint8(postim);
   
   % Save just a part of the slice
   saveim = saveim(cropCoords(1):cropCoords(2), cropCoords(3):cropCoords(4));

   % Resize too the output size
   saveim   = imresize(saveim,outputImageFactorSize); %#ok<PFBNS>   
   imwrite(saveim, allnewpaths{cnt});
   
end %End for loop over movie frames

end % End of makeMovieFromHologramSequence

