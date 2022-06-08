classdef VisionFung2013 < handle & mlaif.AbstractFung2013
	%% ABSTRACTFUNG2013 provides abstractions and reusable implementations of  
    %  Edward K Fung and Richard E Carson.  Cerebral blood flow with [15O]water PET studies using 
    %  an image-derived input function and MR-defined carotid centerlines.  
    %  Phys. Med. Biol. 58 (2013) 1903–1923.  doi:10.1088/0031-9155/58/6/1903

	%  $Revision$
 	%  was created 22-Nov-2021 20:56:03 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlaif/src/+mlaif.
 	%% It was developed on Matlab 9.11.0.1809720 (R2021b) Update 1 for MACI64.  Copyright 2021 John Joowon Lee.

    properties (Dependent)
        dx
        dy
        dz
        Nx
        Ny
        Nz
        petBasename
        petDynamic % contains PET dynamic as ImagingContext2 (LARGE)
        petStatic % contains PET static as ImagingContext2
        tag_idif
    end

	methods 
        
        %% GET
        
        function g = get.dx(this)
            g = this.anatomy.nifti.mmppix(1);
        end
        function g = get.dy(this)
            g = this.anatomy.nifti.mmppix(2);
        end
        function g = get.dz(this)
            g = this.anatomy.nifti.mmppix(3);
        end        
        function g = get.Nx(this)
            g = size(this.anatomy, 1);
        end
        function g = get.Ny(this)
            g = size(this.anatomy, 2);
        end
        function g = get.Nz(this)
            g = size(this.anatomy, 3);
        end
        function g = get.petBasename(this)
            if ~isempty(this.petBasename_)
                g = this.petBasename_;
                return
            end
            if ~isempty(this.petStatic)
                str = this.petStatic.fileprefix;
                re = regexp(str, '\S+_trc-(?<basename>\w+)_proc-\S+', 'names');
                g = re.basename;
                return
            end
            if ~isempty(this.petDynamic)
                str = this.petDynamic.fileprefix;
                re = regexp(str, '\S+_trc-(?<basename>\w+)_proc-\S+', 'names');
                g = re.basename;
                return
            end
            g = '';
            return
        end
        function     set.petBasename(this, s)
            assert(ischar(s))
            this.petBasename_ = s;
        end
        function g = get.petDynamic(this)
            if ~isempty(this.petDynamic_)
                g = copy(this.petDynamic_);
                return
            end
            g = [];
        end
        function     set.petDynamic(this, s)
            assert(isa(s, 'mlfourd.ImagingContext2'))
            this.petDynamic_ = s;
        end
        function g = get.petStatic(this)
            if ~isempty(this.petStatic_)
                g = copy(this.petStatic_);
                return
            end
            g = [];
        end
        function     set.petStatic(this, s)
            assert(isa(s, 'mlfourd.ImagingContext2'))
            this.petStatic_ = s;
        end
        function g = get.tag_idif(this)
            if 0 == this.innerRadius
                g = sprintf('_idif_to%i', this.outerRadius);
                return
            end
            g = sprintf('_idif_%ito%i', this.innerRadius, this.outerRadius);
        end

        %%
		  
 		function this = VisionFung2013(varargin)
 			%% ABSTRACTFUNG2013
            %  @param destinationPath is the path for writing outputs.  Default is Ccir559754Bids.destinationPath.  
            %         Must specify project ID & subject ID.
            %  @param corners from fsleyes NIfTI [ x y z; ... ], [ [RS]; [LS]; [RI]; [LI] ].
            %  @param bbBuffer is the bounding box buffer ~ [x y z] in voxels.
            %  @param iterations ~ 80:130.
            %  @param smoothFactor ~ 0.
            %  @param contractBias is the contraction bias for activecontour():  ~[-1 1], bias > 0 contracting.
            %  @param segmentationOnly is logical.
            %  @param segmentationBlur is scalar.
            %  @param segmentationThresh is scalar.
            %  @param ploton is bool for showing IDIFs.
            %  @param plotqc is bool for showing QC.
            %  @param plotdebug is bool for showing information for debugging.
            %  @param plotclose closes plots after saving them.
            %  @param threshqc for buildCenterline.

            this = this@mlaif.AbstractFung2013(varargin{:});
 			this.bids_ = mlvg.Ccir1211Bids(varargin{:}); 

            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'destinationPath', this.destinationPath, @isfolder)
            addParameter(ip, 'bbBuffer', [0 0 0], @(x) isvector(x) || isa(x, 'containers.Map'))
            addParameter(ip, 'selected_tracers', {''}, @(x) ischar(x) || iscell(x))
            parse(ip, varargin{:})
            ipr = ip.Results;
            if isa(ipr.bbBuffer, 'containers.Map')
                ipr.bbBuffer = ipr.bbBuffer(this.bids.subjectFolder);
            end  
            this.bbBuffer = ipr.bbBuffer;

            % gather requirements
            this.buildTimings();
        end

        function this = buildTimings(this)
            %% builds taus, times, timesMid.

            this.taus = containers.Map;
            this.taus('OC') = this.bids.registry.consoleTaus('OC');
            this.taus('OO') = this.bids.registry.consoleTaus('OO');
            this.taus('HO') = this.bids.registry.consoleTaus('HO');
            this.taus('FDG') = this.bids.registry.consoleTaus('FDG');

            this.times = containers.Map;
            for key = this.taus.keys
                 this.times(key{1}) = [0 cumsum(this.taus(key{1}))];
            end

            this.timesMid = containers.Map;
            for key = this.taus.keys
                taus_ = this.taus(key{1});
                times_ = this.times(key{1});
                this.timesMid(key{1}) = times_(1:length(taus_)) + taus_/2;
            end
        end
    end

    %% PROTECTED
    
    properties (Access = protected)
        petBasename_
        petDynamic_
        petStatic_
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
end

