classdef ArterialInputFunction < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %% ARTERIALINPUTFUNCTION
    %  
    %  Created 29-Mar-2022 22:13:52 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlaif/src/+mlaif.
    %  Developed on Matlab 9.12.0.1884302 (R2022a) for MACI64.  Copyright 2022 John J. Lee.
    
    methods (Static)
    end

    properties (Dependent)
        anatomy % mlfourd.ImagingContext2(T1w|tof) aids visualization of centerline
        anatomy1 % aids visualization of centerline in bounding-box enlarged
        arterialAnatClass % text for logging
        cache_file % .mat
        coord1
        coord2
        doriginof % origin[expanded bb] - origin[bb]
        do_centerline_imdilate % logical for calling imdilate() on centerline
        ic_centerline % imaging context representing centerline in bounding-box
        ic_centerline1 % imaging context representing centerline in bounding-box enlarged
        imdilate_scale % voxels used with imdilate_strel
        imdilate_scale_mm % mm internally converted to voxels and used with imdilate_strel
        imdilate_strel % default is 'sphere'
        mmppix % of anatomy
        pc_centerline % pointCloud representing centerline in bounding-box
        pc_centerline1 % pointCloud representing centerline in bounding-box enlarged
        pet_dyn % imaging context
        pet_dyn_on_anatomy % imaging context
        pet_dyn_on_anatomy1 % imaging context in bounding-box enlarged
        pet_earlyt_on_anatomy % imaging context
        pet_earlyt_on_anatomy1 % imaging context in bounding-box enlarged
        product % best achievable IDIF for requested anatomy and pet_dyn
        reregistration % mlvg.Reregistration
        suffix % char array
        use_cache % logical
    end

    methods

        %% GET, SET

        function g = get.anatomy(this)
            g = copy(this.fung2013_.anatomy);
        end
        function g = get.anatomy1(this)
            if ~isempty(this.anatomy1_)
                g = copy(this.anatomy1_);
                return
            end

            bb = mlaif.ArterialBoundingBox(this.anatomy, this.coord1, this.coord2);
            bb = bb.enlarge_box_base();
            this.anatomy1_ = bb.extract();
            this.anatomy1_.fileprefix = ...
                strcat(this.anatomy1_.fileprefix, '_', clientname(true, 2));
            g = this.anatomy1_;
        end
        function g = get.arterialAnatClass(this)
            g = this.fung2013_.arterialAnatClass;
        end
        function g = get.cache_file(this)
            g = fullfile(this.pet_dyn.filepath, ...
                strcat(clientname(true, 2), this.suffix, '.mat'));
        end
        function g = get.coord1(this)
            g = this.fung2013_.coord1;
        end
        function g = get.coord2(this)
            g = this.fung2013_.coord2;
        end
        function g = get.doriginof(this)
            if ~isempty(this.doriginof_)
                g = this.doriginof_;
                return
            end

            bb = mlaif.ArterialBoundingBox(this.pet_earlyt_on_anatomy, this.coord1, this.coord2);
            originof__ = bb.originof;
            bb.enlarge_box_base();
            this.doriginof_ = bb.originof - originof__;
            g = this.doriginof_;
        end
        function g = get.do_centerline_imdilate(this)
            if ~isempty(this.imdilate_scale) && this.imdilate_scale > 0
                g = true;
                return
            end
            g = false;
        end
        function g = get.ic_centerline(this)
            g = copy(this.ic_centerline_);
        end
        function g = get.ic_centerline1(this)
            g = copy(this.ic_centerline1_);
        end
        function g = get.imdilate_scale(this)
            g = round(this.imdilate_scale_mm / min(this.mmppix));
            if g < 1
                g = [];
            end
        end
        function g = get.imdilate_scale_mm(this)
            g = this.fung2013_.imdilate_scale_mm;
        end
        function g = get.imdilate_strel(this)
            g = this.fung2013_.imdilate_strel;
        end
        function g = get.mmppix(this)
            g = this.anatomy.imagingFormat.mmppix;
        end
        function g = get.pc_centerline(this)
            g = this.pc_centerline_;
        end
        function g = get.pc_centerline1(this)
            if ~isempty(this.pc_centerline1_)
                g = this.pc_centerline1_;
                return
            end

            locs = this.pc_centerline.Location - this.doriginof.*this.mmppix;
            this.pc_centerline1_ = pointCloud(locs);
            g = this.pc_centerline1_;
        end
        function g = get.pet_dyn(this)
            g = copy(this.pet_dyn_);
        end
        function g = get.pet_dyn_on_anatomy(this)
            if ~isempty(this.pet_dyn_on_anatomy_)
                g = copy(this.pet_dyn_on_anatomy_);
                return
            end

            this.pet_dyn_on_anatomy_ = this.fung2013_.pet_dyn_on_anatomy( ...
                this.pet_dyn, this.pet_earlyt_on_anatomy_flirt_);
            g = this.pet_dyn_on_anatomy_;
        end
        function g = get.pet_dyn_on_anatomy1(this)
            if ~isempty(this.pet_dyn_on_anatomy1_)
                g = copy(this.pet_dyn_on_anatomy1_);
                return;
            end

            bb = mlaif.ArterialBoundingBox(this.pet_dyn_on_anatomy, this.coord1, this.coord2);
            bb = bb.enlarge_box_base();
            this.pet_dyn_on_anatomy1_ = bb.extract();
            this.pet_dyn_on_anatomy1_.fileprefix = ...
                strcat(this.pet_dyn_on_anatomy1_.fileprefix, '_', clientname(true, 2));
            g = this.pet_dyn_on_anatomy1_;
        end
        function g = get.pet_earlyt_on_anatomy(this)
            if ~isempty(this.pet_earlyt_on_anatomy_)
                g = copy(this.pet_earlyt_on_anatomy_);
                return
            end

            pet_earlyt_ = this.early_times_averaged(this.pet_dyn);
            [this.pet_earlyt_on_anatomy_,this.pet_earlyt_on_anatomy_flirt_] = this.fung2013_.pet_static_on_anatomy(pet_earlyt_);
            g = this.pet_earlyt_on_anatomy_;
        end
        function g = get.pet_earlyt_on_anatomy1(this)
            if ~isempty(this.pet_earlyt_on_anatomy1_)
                g = copy(this.pet_earlyt_on_anatomy1_);
                return
            end

            bb = mlaif.ArterialBoundingBox(this.pet_earlyt_on_anatomy, this.coord1, this.coord2);
            bb = bb.enlarge_box_base();
            this.pet_earlyt_on_anatomy1_ = bb.extract();
            this.pet_earlyt_on_anatomy1_.fileprefix = ...
                strcat(this.pet_earlyt_on_anatomy1_.fileprefix, '_', clientname(true, 2));
            g = this.pet_earlyt_on_anatomy1_;
        end
        function g = get.product(this)
            g = copy(this.product_);
        end
        function g = get.reregistration(this)
            g = copy(this.reregistration_);
        end
        function g = get.suffix(this)
            c1 = this.fung2013_.str_coord1;
            c2 = this.fung2013_.str_coord2;
            try
                g = sprintf('_%s-%s-%sto%s', this.pet_dyn_on_anatomy.fileprefix, this.arterialAnatClass, c1, c2);
            catch
                g = sprintf('_%s-%s-%sto%s', this.pet_dyn.fileprefix, this.arterialAnatClass, c1, c2);
            end
        end
        function g = get.use_cache(this)
            g = this.use_cache_;
        end

        %%

        function this = ArterialInputFunction(varargin)
            %% ARTERIALINPUTFUNCTION 
            %  Args:
            %      fung2013 (mlaif.Fung2013): director of Fung2013 implementation, reference.
            %      segmentation (mlfourd.ImagingContext2): product of mlaif.ArterialSegmentation.
            %      ic_centerline (mlfourd.ImagingContext2): representing centerline in bounding-box.
            %      imdilate_scale_mm (scalar): mm internally converted to voxels and used with imdilate_strel;
            %                                  0|[] disables imdilate sample_input_function*().
            %      imdilate_strel (text): e.g., 'sphere, 'cube', 'cuboid'; default is 'sphere'
            %      pc_centerline (pointCloud): representating centerline in bounding-box.
            %      pet_dyn (mlfourd.ImagingContext2): dynamic in native scanner space.
            %      use_cache (logical): default is true.
            
            ip = inputParser;
            addParameter(ip, "fung2013", [], @(x) isa(x, 'mlaif.Fung2013'));
            addParameter(ip, "segmentation", [], @(x) isa(x, 'mlfourd.ImagingContext2') || isempty(x));
            addParameter(ip, "imdilate_scale_mm", 0.25, @isnumeric);
            addParameter(ip, "imdilate_strel", 'sphere', @istext);
            addParameter(ip, "ic_centerline", [], @(x) isa(x, 'mlfourd.ImagingContext2') || isempty(x));
            addParameter(ip, "pc_centerline", [], @(x) isa(x, 'pointCloud') || isempty(x));
            addParameter(ip, "pet_dyn", [], @(x) isa(x, 'mlfourd.ImagingContext2'));
            addParameter(ip, "use_cache", true, @islogical);
            parse(ip, varargin{:});
            ipr = ip.Results;
            this.fung2013_ = ipr.fung2013;
            this.segmentation_ = ipr.segmentation;
            this.segmentation_.relocateDerivatives();
            this.ic_centerline_ = ipr.ic_centerline;
            this.ic_centerline_.relocateDerivatives();
            this.fung2013_.imdilate_scale_mm = ipr.imdilate_scale_mm;
            this.fung2013_.imdilate_strel = ipr.imdilate_strel;
            this.pc_centerline_ = ipr.pc_centerline;
            this.pet_dyn_ = ipr.pet_dyn;
            this.pet_dyn_.relocateDerivatives()
            this.use_cache_ = ipr.use_cache;

            if this.use_cache_ && isfile(this.cache_file)
                this.load();
            end
        end

        function this = build_input_function(this)
            %% Builds an input function from centerline samples of dynamic PET
            %  Returns:
            %      this.pet_dyn_on_anatomy: is the full dynamic pet registered using early times.
            %      this.product: best achievable IDIF for requested anatomy and pet_dyn.

            if ~isempty(this.product_)
                return
            end 

            if ~this.fung2013_.needs_reregistration                
                % sample input function from dynamic pet, tight bounding box
                aif_ = this.sample_input_function();
            else
                % re-register ic_centerline1 & pc_centerline1 to pet early times as needed
                this = this.reregister_centerline();    
                % sample input function from dynamic pet, expanded bounding box
                aif_ = this.sample_input_function1();
            end

            % assemble product; save cache_file
            this.product_ = mlfourd.ImagingContext2(aif_);
            save(this);
            %plot(this);
            save(this.pet_earlyt_on_anatomy); % for exploratory QC
        end
        function ic = early_times_averaged(this, ic_dyn)
            if this.is_carbon_monoxide(ic_dyn)
                ic = ic_dyn.timeAveraged();
                ic = ic.blurred(2.5); % 2.5 mm blurring specified by Fung & Carson
                return
            end

            bin = this.msktgen(ic_dyn);
            bin = bin.binarized();
            ic_avgxyz = ic_dyn.volumeAveraged(~bin);
            ic_max = dipmax(ic_avgxyz);
            img = ic_avgxyz.imagingFormat.img;
            [~,it10] = max(img > 0.1*ic_max);
            [~,it95] = max(img > 0.95*ic_max);    
            dit = abs(it95-it10);
            ic = ic_dyn.timeAveraged(it10:it95+dit);
            ic = ic.blurred(2.5);
        end
        function tf = is_carbon_monoxide(~, ic)
            fp = lower(ic.fileprefix);
            if startsWith(fp, 'ocdt') || startsWith(fp, 'codt')
                tf = true;
                return
            end
            if contains(fp, 'trc-oc') || contains(fp, 'trc-co')
                tf = true;
                return
            end
            tf = false;
        end
        function this = load(this)
            cf = load(this.cache_file);
            this.anatomy1_ = cf.this.anatomy1_;
            this.doriginof_ = cf.this.doriginof_;
            this.ic_centerline_ = cf.this.ic_centerline_;
            this.ic_centerline1_ = cf.this.ic_centerline1_;
            this.pc_centerline_ = cf.this.pc_centerline_;
            this.pc_centerline1_ = cf.this.pc_centerline1_;
            this.pet_dyn_ = cf.this.pet_dyn_;
            this.pet_dyn_on_anatomy_ = cf.this.pet_dyn_on_anatomy_;
            this.pet_dyn_on_anatomy1_ = cf.this.pet_dyn_on_anatomy1_;
            this.pet_earlyt_on_anatomy_ = cf.this.pet_earlyt_on_anatomy_;
            this.pet_earlyt_on_anatomy_flirt_ = cf.this.pet_earlyt_on_anatomy_flirt_;
            this.pet_earlyt_on_anatomy1_ = cf.this.pet_earlyt_on_anatomy1_;
            this.product_ = cf.this.product_;
            this.reregistration_ = cf.this.reregistration_;
            this.segmentation_ = cf.this.segmentation_;
        end        
        function ic = msktgen(~, ic_dyn)
            ic_avgt = ic_dyn.timeAveraged();
            ic_avgt.save();
            fn = mlfsl.Flirt.msktgen(ic_avgt.fqfn);
            deleteExisting(ic_avgt);
            ic = mlfourd.ImagingContext2(fn);
        end
        function h = plot(this, varargin)
            h = figure;
            plot(this.product, varargin{:}, 'LineWidth', 2);
            title(clientname(true, 2), 'Interpreter', 'none');
            xlabel('frame', 'FontSize', 14);
            ylabel('specific activity (Bq/mL)', 'FontSize', 14);
            this.savefig(h);
        end
        function h = plott(this, timesMid, plotargs)
            arguments
                this mlaif.ArterialInputFunction
                timesMid double = []
                plotargs = {}
            end
            if isempty(timesMid)
                h = plot(this);
                return
            end

            h = figure;
            plot(timesMid, this.product, plotargs{:}, 'LineWidth', 2);
            title(clientname(true, 2), 'Interpreter', 'none');
            xlabel('frame', 'FontSize', 14);
            ylabel('specific activity (Bq/mL)', 'FontSize', 14);
            this.savefig(h);
        end
        function this = reregister_centerline(this)
            %% Re-register centerline & pc_centerline to pet early times on anatomy as needed;
            %  apply bounding boxes.  
            %  Returns:
            %      this.reregistration_:  contains reregistered properties
            %          ic_centerline1, pc_centerline1, tform1, rewards1, and additional details.

            this.reregister_pcshow('tag', 'pre');

            % pointCloud of bb of pet
            thresh = dipmax(this.pet_earlyt_on_anatomy1)/2;
            pc_pet_static1 = pointCloud(this.pet_earlyt_on_anatomy1, 'thresh', thresh); 

            % re-register
            bb = mlaif.ArterialBoundingBox(this.pet_earlyt_on_anatomy, this.coord1, this.coord2);
            bb = bb.enlarge_box_base();
            this.reregistration_ = mlvg.Reregistration(this.pet_earlyt_on_anatomy1, 'bounding_box', bb);
            tform = rigid3d(eye(4));
            this.reregistration_.pcregistermax(tform, this.pc_centerline1, pc_pet_static1);
            this.pc_centerline1_ = this.reregistration_.pc_centerline1; % altered by reregistration
            this.ic_centerline1_ = this.reregistration_.ic_centerline1; % altered by reregistration

            this.reregister_pcshow('tag', 'post');
        end
        function h = reregister_pcshow(this, varargin)
            if isempty(this.reregistration_)
                markers = 'm';
                note = ' before reregistration';
            else
                markers = 'c';
                note = ' after reregistration';
            end

            % pointCloud of pet, expanded bb
            pet_static1__ = this.pet_earlyt_on_anatomy1;
            pc_pet_static1 = pointCloud(pet_static1__, 'thresh', dipmax(pet_static1__)/2); 

            % pointCloud of anatomy, expanded bb
            scale = dipmax(this.pet_earlyt_on_anatomy1)/dipmax(this.anatomy1)/2;
            anatomy1__ = this.anatomy1 * scale;
            pc_anatomy1 = pointCloud(anatomy1__, 'thresh', dipmax(anatomy1__)/2);    

            % debug
            h = figure('Position', [0 0 1200 1200]);
            xlabel('x/mm');
            ylabel('y/mm');
            hold on
            pcshow(pc_anatomy1);
            pcshow(this.pc_centerline1.Location, markers, 'MarkerSize', 5);
            pcshow(pc_pet_static1);
            hold off
            title(strcat(clientname(true, 2), this.suffix, note), 'Interpreter', 'none');
            this.savefig(h, varargin{:});
        end
        function aif = sample_input_function(this)
            %% Sample input function from dynamic pet on anatomy;
            %  apply tight bounding boxes.
            %  Returns:
            %      aif: vector as mlfourd.ImagingContext2

            % imaging context of centerline, confined to likely interior of artery
            bb = mlaif.ArterialBoundingBox(this.pet_earlyt_on_anatomy, this.coord1, this.coord2);
            e = bb.extract();
            threshed = e.thresh(dipmax(e)/2);
            ic_cl = this.ic_centerline;
            ic_cl = ic_cl .* binarized(threshed);

            % dilate imaging context of centerline
            if this.do_centerline_imdilate
                ic_cl = ic_cl.imdilate_bin(strel(this.imdilate_strel, this.imdilate_scale));
            end

            % sample pet
            bb = mlaif.ArterialBoundingBox(this.pet_dyn_on_anatomy, this.coord1, this.coord2);
            e = bb.extract();
            %ic_cl.ensureSingle();
            aif = e.volumeAveraged(ic_cl);
            aif.fileprefix = strcat(clientname(true, 2), this.suffix);
        end
        function aif = sample_input_function1(this)
            %% Sample input function from dynamic pet on anatomy;
            %  apply expanded bounding boxes.
            %  Returns:
            %      aif: vector as mlfourd.ImagingContext2

            % imaging context of centerline, confined to likely interior of artery
            ic_cl1 = this.ic_centerline1;
            threshed = this.pet_earlyt_on_anatomy1.thresh( ...
                dipmax(this.pet_earlyt_on_anatomy1)/2);
            ic_cl1 = ic_cl1 .* binarized(threshed);

            % dilate imaging context of centerline
            if this.do_centerline_imdilate
                ic_cl1 = ic_cl1.imdilate_bin(strel(this.imdilate_strel, this.imdilate_scale));
            end

            % sample pet
            ic_cl1.ensureSingle();
            aif = this.pet_dyn_on_anatomy1.volumeAveraged(ic_cl1);
            aif.fileprefix = strcat(clientname(true, 2), this.suffix);
        end
        function save(this)
            save(this.cache_file, 'this');
        end
        function savefig(this, h, varargin)
            ip = inputParser;
            addParameter(ip, 'tag', '', @istext)
            parse(ip, varargin{:});
            ipr = ip.Results;
            if ~isempty(ipr.tag)
                ipr.tag = strcat('-', ipr.tag);
            end

            fqfp = fullfile(fileparts(this.cache_file), clientname(true, 3));
            saveas(h, sprintf('%s%s%s.fig', fqfp, this.suffix, ipr.tag));
            saveas(h, sprintf('%s%s%s.png', fqfp, this.suffix, ipr.tag));
            close(h);
        end
    end
    
    %% PROTECTED

    properties (Access = protected)
        anatomy1_
        doriginof_
        fung2013_
        ic_centerline_
        ic_centerline1_
        pc_centerline_
        pc_centerline1_
        pet_dyn_
        pet_dyn_on_anatomy_
        pet_dyn_on_anatomy1_
        pet_earlyt_on_anatomy_
        pet_earlyt_on_anatomy_flirt_
        pet_earlyt_on_anatomy1_
        product_
        reregistration_
        segmentation_
        use_cache_
    end

    methods (Access = protected)
        function that = copyElement(this)
            %%  See also web(fullfile(docroot, 'matlab/ref/matlab.mixin.copyable-class.html'))
            
            that = copyElement@matlab.mixin.Copyable(this);
            that.anatomy1_ = copy(this.anatomy1_);
            that.ic_centerline_ = copy(this.ic_centerline_);
            that.ic_centerline1_ = copy(this.ic_centerline1_);
            that.pet_dyn_ = copy(this.pet_dyn_);
            that.pet_dyn_on_anatomy_ = copy(this.pet_dyn_on_anatomy_);
            that.pet_dyn_on_anatomy1_ = copy(this.pet_dyn_on_anatomy1_);
            that.pet_earlyt_on_anatomy_ = copy(this.pet_earlyt_on_anatomy_);
            that.pet_earlyt_on_anatomy_flirt_ = copy(this.pet_earlyt_on_anatomy_flirt_);
            that.pet_earlyt_on_anatomy1_ = copy(this.pet_earlyt_on_anatomy1_);
            that.product_ = copy(this.product_);
            that.reregistration_ = copy(this.reregistration_);
            that.segmentation_ = copy(this.segmentation_);
        end
    end    
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
