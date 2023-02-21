classdef Fung2013 < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %% 
    %  
    %  Created 25-Mar-2022 15:08:41 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlaif/src/+mlaif.
    %  Developed on Matlab 9.12.0.1884302 (R2022a) for MACI64.  Copyright 2022 John J. Lee.
    
    methods (Static)
        function globbed = createFromCoords(varargin)
            %%
            %  Args:
            %      bids (mlpipeline.IBids): e.g., mlvg.Ccir1211Bids with prop pet_dyn_toglob.
            %      coords (struct): e.g., 
            %      coords = struct( ...
            %               't1w', struct( ...
            %                      'left',  struct( ...
            %                               'coord1', [139 153 14], ...
            %                               'coord2', [157 141 117]), ...
            %                      'right', struct( ...
            %                               'coord1', [76 154 14], ...
            %                               'coord2', [62 141 117])), ...
            %               'tof', struct( ...
            %                      'left',  struct( ...
            %                               'coord1', [500 386 37], ...
            %                               'coord2', [388 520 90]), ...
            %                      'right', struct( ...
            %                               'coord1', [228 392 37], ...
            %                               'coord2', [358 473 90])));
            %  Returns:
            %      globbed: cell array of globbed PET dynamic.

            import mlaif.Fung2013;

            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'bids', [], @(x) isa(x, 'mlpipeline.IBids'));
            addParameter(ip, 'coords', [], @isstruct);
            parse(ip, varargin{:});
            ipr = ip.Results;

            globbed = globT(ipr.bids.pet_dyn_toglob);
            for t = {'tof', 't1w'}
                t_ = t{1};
                for h = {'left', 'right'}
                    h_ = h{1};
                    for ig = 1:length(globbed)
                        try
                            f = Fung2013.(strcat('createFor', upper(t_(1)),t_(2:end)))( ...
                                'coord1', ipr.coords.(t_).(h_).coord1, ...
                                'coord2', ipr.coords.(t_).(h_).coord2, ...
                                'use_cache', true, varargin{:});
                            call(f, 'pet_dyn', globbed{ig});
                        catch ME
                            handwarning(ME);
                        end
                    end
                end
            end
        end
        function this = createForT1w(varargin)
            %% Args:
            %      bids (mlpipeline.IBids)
            %      verbose (scalar): 0 is silent; default 1 provides QC; 2+ aids debugging.

            aa = mlaif.ArterialAnatomy.createForT1w(varargin{:});
            this = mlaif.Fung2013('arterial_anatomy', aa, varargin{:});
        end
        function this = createForTof(varargin)
            %% Args:
            %      bids (mlpipeline.IBids)
            %      verbose (scalar): 0 is silent; default 1 provides QC; 2+ aids debugging.

            aa = mlaif.ArterialAnatomy.createForTof(varargin{:});
            this = mlaif.Fung2013('arterial_anatomy', aa, varargin{:});
        end
        function this = createForTofOnT1w(varargin)
            %% Args:
            %      bids (mlpipeline.IBids)
            %      verbose (scalar): 0 is silent; default 1 provides QC; 2+ aids debugging.

            aa = mlaif.ArterialAnatomy.createForTofOnT1w(varargin{:});
            this = mlaif.Fung2013('arterial_anatomy', aa, varargin{:});
        end
        function viewAnat(varargin)
            %% 
            % Args:
            %      bids (mlpipeline.IBids)

            ip = inputParser;
            addRequired(ip, 'bids', [], @(x) isa(x, 'mlpipeline.IBids'));
            parse(ip, varargin{:});
            ipr = ip.Results;

            ipr.bids.t1w_ic.view();
            ipr.bids.tof_ic.view();
        end
    end

    properties (Dependent)
        anatomy
        arterialAnatClass % text for logging
        bids
        coord1            % ijk of bounding box
        coord2            % ijk of bounding box
        exclusion_blur    % numerical, depends on ECIC | ICIC
        exclusion_init    % mlfourd.ImagingContext2
        exclusion_list    % cell
        exclusion_thresh  % numerical
        imdilate_scale_mm % mm internaly converted to voxels and used with imdilate_strel
        imdilate_strel    % default is 'sphere'
        inclusion_blur    % numerical, depends on ECIC | ICIC
        inclusion_init    % mlfourd.ImagingContext2
        inclusion_thresh  % numerical
        needs_reregistration % logical
        pc_threshp
        product
        pet_dyn           % mlfourd.ImagingContext2
        str_coord1        % as nice char array
        str_coord2        % as nice char array
        verbose           % 0 is silent; default 1 provides QC; 2+ aids debugging
        wmparc

        %% for B-splines in mlvg.Hunyadi2021

        k 
        t
        N_centerline_samples
    end

    methods

        %% GET

        function g = get.anatomy(this)
            g = this.arterial_anatomy_.anatomy;
        end
        function g = get.arterialAnatClass(this)
            g = class(this.arterial_anatomy_);
            g = strsplit(g, '.');
            g = g{end};
        end
        function g = get.bids(this)
            g = this.arterial_anatomy_.bids;
        end
        function g = get.coord1(this)
            g = this.coord1_;
        end
        function g = get.coord2(this)
            g = this.coord2_;
        end
        function g = get.exclusion_blur(this)
            g = this.arterial_anatomy_.exclusion_blur;
        end
        function g = get.exclusion_init(this)
            g = this.arterial_anatomy_.exclusion_init;
        end
        function g = get.exclusion_list(this)
            g = this.arterial_anatomy_.exclusion_list;
        end
        function g = get.exclusion_thresh(this)
            g = this.arterial_anatomy_.exclusion_thresh;
        end
        function g = get.imdilate_scale_mm(this)
            g = this.arterial_anatomy_.imdilate_scale_mm;
        end
        function     set.imdilate_scale_mm(this, s)
            this.arterial_anatomy_.imdilate_scale_mm = s;
        end
        function g = get.imdilate_strel(this)
            g = this.arterial_anatomy_.imdilate_strel;
        end
        function     set.imdilate_strel(this, s)
            this.arterial_anatomy_.imdilate_strel = s;
        end
        function g = get.inclusion_blur(this)
            g = this.arterial_anatomy_.inclusion_blur;
        end
        function g = get.inclusion_init(this)
            g = this.arterial_anatomy_.inclusion_init;
        end
        function g = get.inclusion_thresh(this)
            g = this.arterial_anatomy_.inclusion_thresh;
        end
        function g = get.needs_reregistration(this)
            g = ~isa(this.arterial_anatomy_, 'mlaif.ICIC');
        end
        function g = get.pc_threshp(this)
            g = this.pc_threshp_;
        end
        function g = get.pet_dyn(this)
            g = copy(this.pet_dyn_);
        end
        function     set.pet_dyn(this, s)
            this.pet_dyn_ = mlfourd.ImagingContext2(s);
        end
        function g = get.product(this)
            g = this.arterial_input_function_.product;
        end
        function g = get.str_coord1(this)
            g = strrep(strrep(strrep(mat2str(this.coord1), ' ', ','), '[', ''), ']', '');
        end
        function g = get.str_coord2(this)
            g = strrep(strrep(strrep(mat2str(this.coord2), ' ', ','), '[', ''), ']', '');
        end
        function g = get.verbose(this)
            g = this.verbose_;
        end
        function g = get.wmparc(this)
            g = this.arterial_anatomy_.wmparc;
        end

        function g = get.k(this)
            if isa(this.arterial_anatomy_, 'mlaif.ECIC')
                g = 4; % cubic b-splines
                return
            end
            if isa(this.arterial_anatomy_, 'mlaif.ICIC')
                g = 5; % quartic b-splines
                return
            end
        end
        function g = get.t(this)
            if isa(this.arterial_anatomy_, 'mlaif.ECIC')
                g = [0 0 0 0 0.2 0.4 0.6 0.8 1 1 1 1];
                return
            end
            if isa(this.arterial_anatomy_, 'mlaif.ICIC')
                g = [0 0 0 0 0 0.5 0.7 1 1 1 1 1];
                return
            end
        end
        function g = get.N_centerline_samples(this)
            if isa(this.arterial_anatomy_, 'mlaif.ECIC')
                g = abs(this.coord2(3) - this.coord1(3)) + 1; % sample voxel barycenters
                return
            end
            if isa(this.arterial_anatomy_, 'mlaif.ICIC')
                g = ceil(8*norm(this.coord2 - this.coord1));
                return
            end
        end

        %%

        function call(this, varargin)
            %  Args:
            %      coord1 (vector):  updates values from ctor.
            %      coord2 (vector):  updates values from ctor.
            %      use_cache (logical):  reuse cache file for ArterialSegmentation and ArterialCenterline, but
            %                            not for ArterialInputFunction.

            ip = inputParser;
            addParameter(ip, 'coord1', this.coord1, @isvector);
            addParameter(ip, 'coord2', this.coord2, @isvector);
            addParameter(ip, 'pet_dyn', this.pet_dyn, @(x) ~isempty(x));
            addParameter(ip, 'use_cache', true, @islogical);
            parse(ip, varargin{:})
            ipr = ip.Results;
            this.coord1_ = ipr.coord1;
            this.coord2_ = ipr.coord2;
            this.pet_dyn_ = mlfourd.ImagingContext2(ipr.pet_dyn);

            this.arterial_segmentation_ = mlaif.ArterialSegmentation( ...
                'fung2013', this, ...
                'use_cache', ipr.use_cache);
            this.arterial_segmentation_.build_segmentation();

            this.arterial_centerline_ = mlaif.ArterialCenterline( ...
                'fung2013', this, ...
                'segmentation', this.arterial_segmentation_.product, ...
                'use_cache', ipr.use_cache);
            this.arterial_centerline_.build_centerline();

            this.arterial_input_function_ = mlaif.ArterialInputFunction( ...
                'fung2013', this, ...
                'segmentation', this.arterial_segmentation_.product, ...
                'ic_centerline', this.arterial_centerline_.product, ...
                'pc_centerline', this.arterial_centerline_.pointCloud_for_centerline(), ...
                'pet_dyn', this.pet_dyn, ...
                'use_cache', false);            
            this.arterial_input_function_.build_input_function();

            aif = this.arterial_input_function_;
            save(aif.product);
            save(aif.pet_earlyt_on_anatomy);
            save(aif.ic_centerline);
        end
        function c = complement_coord(this, c, idx)
            %  Args:
            %      c (vector): ijk for voxel.
            %      idx (scalar): coord element to complement.

            assert(isvector(c))
            assert(isscalar(idx))
            sz = size(this.anatomy);
            c(idx) = sz(idx) - c(idx) + 1;           
        end
        function [ic,f] = pet_static_on_anatomy(this, pet)
            [ic,f] = this.arterial_anatomy_.pet_static_on_anatomy(pet);
        end
        function ic = pet_dyn_on_anatomy(this, pet_dyn, f)
            ic = this.arterial_anatomy_.pet_dyn_on_anatomy(pet_dyn, f);
        end
    end

    %% PROTECTED

    properties (Access = protected)
        arterial_anatomy_
        arterial_segmentation_
        arterial_centerline_
        arterial_input_function_
        coord1_
        coord2_
        pc_threshp_
        pet_dyn_
        verbose_
    end

    methods (Access = protected)
        function this = Fung2013(varargin)
            %% FUNG2013 
            %  Args:
            %      arterial_anatomy (mlaif.ArterialAnatomy): factory|strategy picks ECIC|ICIC|others.
            %      coord1 (vector): for bounding box, arranged so that norm(coord1) < norm(coord2).
            %      coord2 (vector): for bounding box, arranged so that norm(coord1) < norm(coord2).
            %      pc_threshp (scalar): default 40% of pointCloud intensities.  
            %      pet_dyn (any): provides dynamic PET for Fung's 2013 method.
            %      verbose (scalar): 0 is silent; 1 provides QC; 2+ aids debugging.
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, "arterial_anatomy", [], @(x) isa(x, 'mlaif.ArterialAnatomy'));
            addParameter(ip, "coord1", [], @isvector);
            addParameter(ip, "coord2", [], @isvector);
            addParameter(ip, "pc_threshp", 35, @isscalar);
            addParameter(ip, "pet_dyn", []);
            addParameter(ip, "verbose", 1, @isscalar);
            parse(ip, varargin{:});
            ipr = ip.Results;
            
            this.arterial_anatomy_ = ipr.arterial_anatomy;
            if norm(ipr.coord1) < norm(ipr.coord2)
                this.coord1_ = ipr.coord1;
                this.coord2_ = ipr.coord2;
            else
                this.coord1_ = ipr.coord2;
                this.coord2_ = ipr.coord1;
            end
            this.pc_threshp_ = ipr.pc_threshp;
            this.pet_dyn_ = mlfourd.ImagingContext2(ipr.pet_dyn);
            this.verbose_ = ipr.verbose;
        end
        function that = copyElement(this)
            %%  See also web(fullfile(docroot, 'matlab/ref/matlab.mixin.copyable-class.html'))
            
            that = copyElement@matlab.mixin.Copyable(this);
            that.arterial_anatomy_ = copy(this.arterial_anatomy_);
            that.arterial_segmentation_ = copy(this.arterial_segmentation_);
            that.arterial_centerline_ = copy(this.arterial_centerline_);
            that.arterial_input_function_ = copy(this.arterial_input_function_);
            that.pet_dyn_ = copy(this.pet_dyn_);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
