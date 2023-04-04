classdef ArterialSegmentation < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %% ARTERIALSEGMENTATION
    %  >> as = mlaif.ArterialSegmentation('fung2013', fung)
    %  >> as.build_segmentation()
    %  >> patch(as)
    %  >> pcshow(as)
    %  
    %  Created 06-Mar-2022 23:35:52 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlaif/src/+mlaif.
    %  Developed on Matlab 9.11.0.1873467 (R2021b) Update 3 for MACI64.  Copyright 2022 John J. Lee.

    properties
        contractBias % open to tweaking
        iterations   %
        smoothFactor %
    end

    properties (Dependent)
        anatomy
        arterialAnatClass % text for logging
        cache_file % .mat
        coord1
        coord2
        exclusion_blur   % numerical, depends on ECIC | ICIC
        exclusion_init   % mlfourd.ImagingContext2
        exclusion_list   % cell
        exclusion_thresh % numerical
        inclusion_blur   % numerical, depends on ECIC | ICIC
        inclusion_init   % mlfourd.ImagingContext2
        inclusion_thresh % numerical
        pc_threshp
        product
        suffix
        use_cache % logical
        wmparc
    end

    methods

        %% GET

        function g = get.anatomy(this)
            g = copy(this.fung2013_.anatomy);
        end
        function g = get.arterialAnatClass(this)
            g = this.fung2013_.arterialAnatClass;
        end
        function g = get.cache_file(this)
            g = fullfile(this.anatomy.filepath, ...
                strcat(clientname(true, 2), this.suffix, '.mat'));
        end
        function g = get.coord1(this)
            g = this.fung2013_.coord1;
        end
        function g = get.coord2(this)
            g = this.fung2013_.coord2;
        end
        function g = get.exclusion_blur(this)
            g = this.fung2013_.exclusion_blur;
        end
        function g = get.exclusion_init(this)
            g = this.fung2013_.exclusion_init;
        end
        function g = get.exclusion_list(this)
            g = this.fung2013_.exclusion_list;
        end
        function g = get.exclusion_thresh(this)
            g = this.fung2013_.exclusion_thresh;
        end
        function g = get.inclusion_blur(this)
            g = this.fung2013_.inclusion_blur;
        end
        function g = get.inclusion_init(this)
            g = this.fung2013_.inclusion_init;
        end
        function g = get.inclusion_thresh(this)
            g = this.fung2013_.inclusion_thresh;
        end
        function g = get.pc_threshp(this)
            g = this.fung2013_.pc_threshp;
        end
        function g = get.product(this)
            g = copy(this.product_);
        end
        function g = get.suffix(this)
            c1 = this.fung2013_.str_coord1;
            c2 = this.fung2013_.str_coord2;
            g = sprintf('_%s-%sto%s', this.arterialAnatClass, c1, c2);
        end
        function g = get.use_cache(this)
            g = this.use_cache_;
        end
        function g = get.wmparc(this)
            g = this.fung2013_.wmparc;
        end

        %%

        function this = ArterialSegmentation(varargin)
            %  Args:
            %      fung2013 (mlaif.Fung2013): director of Fung2013 implementation, reference
            %      contractBias (scalar)
            %      iterations (scalar)
            %      smoothFactor (scalar)  
            %      use_cache (logical): default is true.   
            
            ip = inputParser;
            addParameter(ip, "fung2013", [], @(x) isa(x, 'mlaif.Fung2013'))
            addParameter(ip, "contractBias", 0.2, @isscalar);
            addParameter(ip, "iterations", 80, @isscalar);
            addParameter(ip, "smoothFactor", 0, @isscalar);
            addParameter(ip, "use_cache", true, @islogical);
            parse(ip, varargin{:});
            ipr = ip.Results;
            
            this.fung2013_ = ipr.fung2013;
            this.contractBias = ipr.contractBias;
            this.iterations = ipr.iterations;
            this.smoothFactor = ipr.smoothFactor;
            this.use_cache_ = ipr.use_cache;

            if this.use_cache_ && isfile(this.cache_file)
                this.load();
            end
        end

        function [e,bb] = build_extracted(this)
            %% Simply selects boundary box from this.anatomy.
            %  Returns:
            %      e (mlfourd.ImagingContext2): size delimited by coord1, coord2
            %      bb (mlaif.ArterialBoundingBox)

            bb = mlaif.ArterialBoundingBox(this.anatomy, this.coord1, this.coord2);
            e = bb.extract();
        end
        function [em,bb] = build_extracted_mask(this)
            %% Gathers exclusions from this.wmparc and inclusions from this.anatomy.
            %  Returns:
            %      em (mlfourd.ImagingContext2): size delimited by coord1, coord2
            %      bb (mlaif.ArterialBoundingBox)
            
            exclusion = this.exclusion_init;
            for ie = 1:length(this.exclusion_list)
                e = this.wmparc.select_roi(this.exclusion_list{ie}); % mask out exclusion_list requested of wmparc
                e = e.blurred(this.exclusion_blur(ie));
                e = e.thresh(this.exclusion_thresh(ie));
                e = e.binarized();
                exclusion = exclusion | e;
            end

            inclusion = this.inclusion_init.blurred(this.inclusion_blur);
            inclusion = inclusion.thresh(this.inclusion_thresh);
            inclusion = inclusion.binarized();

            %assert(all(inclusion.imagingFormat.originator == exclusion.imagingFormat.originator))
            disp(exclusion.imagingFormat.hdr.hist)
            disp(inclusion.imagingFormat.hdr.hist)
            bb = mlaif.ArterialBoundingBox(inclusion .* ~exclusion, this.coord1, this.coord2);
            em = bb.extract();
        end
        function this = build_segmentation(this, varargin)
            %% Segments the arterial structure using activecontour() with the 'Chan-Vese' method.
            %  Args:
            %      iterations (scalar): ~100.
            %      method (text): 'Chan-Vese' | 'edge'.
            %      contractBias (-1 < scalar < 1): ~0.02 for 'Chan-Vese'; ~0.3 for 'edge'.
            %          Tendency of the contour to grow outwards or shrink inwards, specified as the comma-separated 
            %          pair consisting of 'ContractionBias' and a numeric scalar. Positive values bias the contour to 
            %          shrink inwards (contract). Negative values bias the contour to grow outwards (expand). This 
            %          parameter does not guarantee that the contour contracts or expands. It is possible that even 
            %          with a positive value for this parameter, the contour could actually expand. However, by 
            %          specifying a bias, you slow the expansion when compared to an unbiased contour.
            %      smoothFactor (scalar >= 0): ~ 0 for 'Chan-Vese'; ~ 1 for 'edge'.
            %  Returns:
            %      this.product: mlfourd.ImagingContext2 of results from segmentation methods.
            
            ip = inputParser;
            addParameter(ip, 'iterations', this.iterations, @isscalar);
            addParameter(ip, 'method', 'Chan-Vese', @isscalar);
            addParameter(ip, 'contractBias', this.contractBias, @isscalar);
            addParameter(ip, 'smoothFactor', this.smoothFactor, @isscalar);
            parse(ip, varargin{:});
            ipr = ip.Results;

            if ~isempty(this.product_)
                return
            end                

            % call Chan-Vese snakes; iterate
            [e,bb] = this.build_extracted();
            em = this.build_extracted_mask();
            ac = activecontour(e.imagingFormat.img, em.imagingFormat.img,  ipr.iterations, ipr.method, ...
                'ContractionBias', ipr.contractBias, ...
                'SmoothFactor', ipr.smoothFactor);

            % product := segmentation in bounding box
            ic = this.anatomy.zeros;
            ic.fileprefix = strcat(clientname(true, 2), this.suffix);
            nii = ic.imagingFormat;
            nii.img = ac;
            
            %% KLUDGE
            %  N.B. mlfourd.ImagingFormatTool.set.originator, mlfourd_unittest.Test_ImagingContext2
            nii.originator(1:3) = ic.imagingFormat.originator(1:3) - bb.originof + ic.qfac*[2,0,0]./bb.mmppix; 
            this.product_ = mlfourd.ImagingContext2(nii);
            save(this);
            %patch(this);
            %pcshow(this);
        end
        function this = load(this)
            cf = load(this.cache_file);
            this.product_ = cf.this.product_;    
            this.contractBias = cf.this.contractBias;
            this.iterations   = cf.this.iterations;
            this.smoothFactor = cf.this.smoothFactor;
        end
        function ax = pcshow(this, varargin)
            %% Embeds the segmentation product in this.anatomy, then pcshows.

            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'threshp', this.pc_threshp, @isscalar)
            parse(ip, varargin{:});
            ipr = ip.Results;

            if isempty(this.product)
                h = figure;
                ax = pcshow(this.anatomy, 'threshp', ipr.threshp, varargin{:});
                return
            end
            abb = mlaif.ArterialBoundingBox(this.anatomy, this.coord1, this.coord2);
            bb = abb.insert(this.product .* dipmax(this.anatomy));
            h = figure;
            ax = pcshow(bb, 'threshp', ipr.threshp, varargin{:});
            title(this);
            this.savefig(h)
        end
        function p = patch(this, varargin)
            %% Plots segmentation rendered as polygonal regions.

            if isempty(this.product)
                figure;
                p = patch(this.anatomy, varargin{:});
                return
            end
            h = figure;
            p = patch(this.product, varargin{:});
            title(this);
            this.savefig(h)
        end
        function save(this)
            save(this.cache_file, 'this');
        end
        function savefig(this, h)
            fqfp = fullfile(fileparts(this.cache_file), clientname(true, 3));
            saveas(h, sprintf('%s%s.fig', fqfp, this.suffix));
            saveas(h, sprintf('%s%s.png', fqfp, this.suffix));
            close(h);
        end
        function title(this)
            title(sprintf( ...
                "%s:\n" + ...
                "iterations %i, contract bias %g, smooth factor %g,\n" + ...
                "excl. blur %g, excl. list %s, excl. thresh %g,\n" + ...
                "incl. blur %g, incl. fileprefix %s, incl. thresh", ...
                strcat(clientname(true, 3), this.suffix), ...
                this.iterations, this.contractBias, this.smoothFactor, ...
                this.exclusion_blur, cell2str(this.exclusion_list), this.exclusion_thresh, ...
                this.inclusion_blur, this.inclusion_init.fileprefix, this.inclusion_thresh), ...
                'Interpreter', 'none', ...
                'FontSize', 10)
        end
    end

    %% PROTECTED

    properties (Access = protected)
        fung2013_
        product_
        use_cache_
    end

    methods (Access = protected)
        function that = copyElement(this)
            %%  See also web(fullfile(docroot, 'matlab/ref/matlab.mixin.copyable-class.html'))
            
            that = copyElement@matlab.mixin.Copyable(this);
            that.product_ = copy(this.product_);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
