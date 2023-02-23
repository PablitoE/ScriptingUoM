classdef PropagationHandler < handle
   % This class simplifies the propagation of a hologram to a single z
   % position. The working hologram may or not change but the config file
   % remains always the same.
   properties
      hologramFilename = '';
      configPath = '';
      cfg = [];
      hologram =[];
      propagator = [];
      slice = [];
      patch = [];
      patchthr = [];
   end

   methods
      function this = PropagationHandler(cfgPath)
         this.configPath = cfgPath;
         this.cfg = config(this.configPath);
      end

      function updateSlice(this, newHologramFilename, z)
         if ~strcmp(this.hologramFilename,newHologramFilename)
            this.hologramFilename = newHologramFilename;
            this.cfg.current_holo = newHologramFilename;
            % Prepare the hologram
            preimage = img;
            preimage.filter_handle = 'prefilters';
            preimage.config_handle = this.cfg;
            preimage.unregisterListeners;
            this.hologram = preimage.ampFilter(this.cfg.full_holo_path);
            preimage.delete;
            clear preimage;

            this.propagator = Propagator(); % Propagation object
            this.propagator.should_normalize = true;
            this.propagator.config_handle = this.cfg;
            if isempty(gcp('nocreate'))
               if isdeployed,  this.propagator.fft_workers = 1; else
                  num_logical_processors = str2double(getenv('NUMBER_OF_PROCESSORS'));
                  this.propagator.fft_workers = floor(num_logical_processors / 4);  % Use 2 workers for small images or 6 workers for larger images (less processors due to memory reqs)
               end
            else, this.propagator.fft_workers = 2;
            end
            this.propagator.verbose = false;
            this.propagator.preconstruct(this.hologram,z); % With this hologram
            if ~this.cfg.getDynParam('useSAASM'), this.propagator.FPrepped_root;
               if ~this.cfg.getDynParam('useBLASM'), this.propagator.FPrepped_filter; end
            end
            this.propagator.unregisterListeners;
         end

         % Do the reconstruction/propagation
         [this.slice, ~]   = this.propagator.slice(z);
      end

      function rescalePatch(this,rows,cols)
         patch_original = this.slice(rows(1):rows(2),cols(1):cols(2));
         [this.patchthr,this.patch] = particles.rescaleImage(patch_original,this.cfg.getDynParam('rescaleThresh'),this.cfg.closeGapsRad,this.cfg.shouldFillHoles);
      end

      function centroid = patchCentroid(this)
         if ~isempty(this.patchthr)
            s = regionprops(this.patchthr,'centroid');
            centroid = s.Centroid;
         else
            centroid = [];
         end
      end
   end
end