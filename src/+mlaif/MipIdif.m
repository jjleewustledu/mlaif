classdef MipIdif < handle & mlsystem.IHandle
    %% is a builder.  Start with build_all(this).  
    %  
    %  Created 12-Jun-2023 22:59:36 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlaif/src/+mlaif.
    %  Developed on Matlab 9.14.0.2254940 (R2023a) Update 2 for MACI64.  Copyright 2023 John J. Lee.
    
    properties (Dependent)
        centerline_on_pet
        centerline_source
        halflife
        pet_avgt % w/ decay-correction, weighted by taus
        pet_avgt_on_tof
        pet_dyn % read from filesystem to minimize copies of objects in memory, relocated from sourcedata -> derivatives
        pet_mipt % w/ decay-correction
        pet_mipt_on_tof
        product
        t1w
        tof
        tof_on_pet
        tof_on_pet_flirt
        tof_on_t1w
        tracer
        weights_timesMid
    end

    methods %% GET, SET
        function g = get.centerline_on_pet(this)
            if isempty(this.centerline_on_pet_)
                try
                    this.centerline_on_pet_ = mlfourd.ImagingContext2( ...
                        fullfile(this.pet_dyn.filepath, "centerline_on_pet.nii.gz"));
                catch ME
                    handexcept(ME)
                end
            end
            g = copy(this.centerline_on_pet_);
        end
        function g = get.centerline_source(this)
            globbed = glob( ...
                convertStringsToChars( ...
                fullfile(this.centerline_on_pet.filepath, "sub-*_ses-*_trc-*Static.nii.gz")));
            assert(~isempty(globbed), stackstr())
            g = mlfourd.ImagingContext2(globbed{1});
        end
        function g = get.halflife(this)
            rn = this.tracer_kit_.make_radionuclides();
            g = rn.halflife;
            assert(isfinite(g))
        end
        function g = get.pet_avgt(this)
            if ~isempty(this.pet_avgt_)
                g = this.pet_avgt_;
                return
            end

            fqfn = strrep(this.pet_dyn.fqfp, "createNiftiMovingAvgFrames", "createNiftiStatic");
            if isfile(fqfn)
                this.pet_avgt_ = mlfourd.ImagingContext2(fqfn);
                g = this.pet_avgt_;
                return
            end

            fqfn = this.pet_dyn.fqfp+"_avgt.nii.gz";
            if isfile(fqfn)
                this.pet_avgt_ = mlfourd.ImagingContext2(fqfn);
                g = this.pet_avgt_;
                return
            end

            this.pet_avgt_ = this.pet_dyn.timeAveraged(weights=this.weights_timesMid);
            this.pet_avgt_.relocateToDerivativesFolder();
            save(this.pet_avgt_);
            g = this.pet_avgt_;
        end
        function g = get.pet_avgt_on_tof(this)
            if ~isempty(this.pet_avgt_on_tof_)
                g = this.pet_avgt_on_tof_;
                return
            end

            fqfn = this.pet_avgt.fqfp+"_on_tof.nii.gz";
            if isfile(fqfn)
                this.pet_avgt_on_tof_ = mlfourd.ImagingContext2(fqfn);
                g = this.pet_avgt_on_tof_;
                return
            end

            this.build_pet_objects();
            g = this.pet_avgt_on_tof_;
        end
        function g = get.pet_dyn(this)
            g = mlfourd.ImagingContext2(this.pet_dyn_fqfn_);
        end
        function g = get.pet_mipt(this)
            if ~isempty(this.pet_mipt_)
                g = this.pet_mipt_;
                return
            end

            fqfn = this.pet_dyn.fqfp+"_mipt.nii.gz";
            if isfile(fqfn)
                this.pet_mipt_ = mlfourd.ImagingContext2(fqfn);
                g = this.pet_mipt_;
                return
            end

            ifc = this.pet_dyn.imagingFormat;
            weights = this.weights_timesMid;
            for t = 1:size(ifc, 4)
                ifc.img(:,:,:,t) = ifc.img(:,:,:,t)*weights(t);
            end

            ic = mlfourd.ImagingContext2(ifc);
            this.pet_mipt_ = max(ic, [], 4);
            this.pet_mipt_.fileprefix = sprintf("%s_mipt", this.pet_dyn.fileprefix);
            this.pet_mipt_.relocateToDerivativesFolder();
            save(this.pet_mipt_);
            g = this.pet_mipt_;
        end
        function g = get.pet_mipt_on_tof(this)
            if ~isempty(this.pet_mipt_on_tof_)
                g = this.pet_mipt_on_tof_;
                return
            end

            fqfn = this.pet_mipt.fqfp+"_on_tof.nii.gz";
            if isfile(fqfn)
                this.pet_mipt_on_tof_ = mlfourd.ImagingContext2(fqfn);
                g = this.pet_mipt_on_tof_;
                return
            end

            this.build_pet_objects();
            g = this.pet_mipt_on_tof_;
        end
        function g = get.product(this)
            assert(~isempty(this.product_))
            g = copy(this.product_);
        end
        function g = get.t1w(this)
            med = this.bids_kit_.make_bids_med();
            g = med.t1w_ic;
        end
        function g = get.tof(this)
            med = this.bids_kit_.make_bids_med();
            g = med.tof_ic;
        end
        function g = get.tof_on_pet(this)
            if ~isempty(this.tof_on_pet_)
                g = this.tof_on_pet_;
                return
            end

            fn = this.tof.fileprefix+"_on_pet_avgt.nii.gz";
            fqfn = fullfile(this.pet_dyn.filepath, fn);
            if isfile(fqfn)
                this.tof_on_pet_ = mlfourd.ImagingContext2(fqfn);
                g = this.tof_on_pet_;
                return
            end

            this.build_pet_objects();
            g = this.tof_on_pet_;
        end
        function g = get.tof_on_pet_flirt(this)
            if ~isempty(this.tof_on_pet_flirt_)
                g = this.tof_on_pet_flirt_;
                return
            end
            this.build_pet_objects();
            g = this.tof_on_pet_flirt_;
        end
        function g = get.tof_on_t1w(this)
            med = this.bids_kit_.make_bids_med();
            g = med.tof_on_t1w_ic;
        end
        function g = get.tracer(this)
            if ~isempty(this.tracer_)
                g = this.tracer_;
                return
            end

            med = this.bids_kit_.make_bids_med();
            g = med.tracer;
            this.tracer_ = g;
        end
        function g = get.weights_timesMid(this)
            if ~isempty(this.weights_timesMid_)
                g = this.weights_timesMid_;
                return
            end

            j = this.pet_dyn.json_metadata;
            timesMid = j.timesMid;
            timesMid = mlpipeline.ImagingMediator.ensureNumericTimesMid(timesMid);
            timesMid = asrow(timesMid);
            try
                weights = exp(-timesMid/this.halflife);
                weights(~isfinite(weights)) = 0;
                this.weights_timesMid_ = weights/sum(weights); % normalize to unit AUC
                g = this.weights_timesMid_;
                assert(all(isfinite(g)))
            catch ME
                handwarning(ME)
                g = ones(size(timesMid));
            end
        end
    end

    methods
        function cl = back_project(this)
            %% back-project MIPs            
            
            pwd1 = pushd(this.pet_dyn.filepath);

            % read_avw <- $FSLDIR/etc/matlab
            cxL_img = read_avw(fullfile(this.pet_dyn.filepath, "centerline_xl.nii.gz"));
            cxR_img = read_avw(fullfile(this.pet_dyn.filepath, "centerline_xr.nii.gz"));
            cyL_img = read_avw(fullfile(this.pet_dyn.filepath, "centerline_yl.nii.gz"));
            cyR_img = read_avw(fullfile(this.pet_dyn.filepath, "centerline_yr.nii.gz"));
            czL_img = read_avw(fullfile(this.pet_dyn.filepath, "centerline_zl.nii.gz"));
            czR_img = read_avw(fullfile(this.pet_dyn.filepath, "centerline_zr.nii.gz"));
            
            size_tof = size(this.tof);
            img1 = zeros(size_tof);
            img2 = zeros(size_tof);
            img3 = zeros(size_tof);
            img4 = zeros(size_tof);
            img5 = zeros(size_tof);
            img6 = zeros(size_tof);
            
            for ix = 1:size_tof(1)
                img1(ix,:,:) = cxL_img; end
            for iy = 1:size_tof(2)
                img2(:,iy,:) = cyL_img; end
            for iz = 1:size_tof(3)
                img3(:,:,iz) = czL_img; end
            for ix = 1:size_tof(1)
                img4(ix,:,:) = cxR_img; end
            for iy = 1:size_tof(2)
                img5(:,iy,:) = cyR_img; end
            for iz = 1:size_tof(3)
                img6(:,:,iz) = czR_img; end
            
            img1 = imdilate(img1, strel("sphere", 1));
            img2 = imdilate(img2, strel("sphere", 1));
            img3 = imdilate(img3, strel("sphere", 1));
            img4 = imdilate(img4, strel("sphere", 1));
            img5 = imdilate(img5, strel("sphere", 1));
            img6 = imdilate(img6, strel("sphere", 1));
            
            imgL = img1 & img2 & img3;
            imgR = img4 & img5 & img6;
            
            %tof1 = this.tof.thresh(150).binarized();
            cl = copy(this.tof); 
            cl = cl.nifti;
            cl.img = imgL + imgR;
            cl.fileprefix = 'centerline_on_tof';
            cl = mlfourd.ImagingContext2(cl);
            cl.save()
            %cl.pcshow()  

            % transform centerline from tof native to pet native
            this.centerline_on_pet_ = mlfourd.ImagingContext2( ...
                fullfile(this.pet_dyn.filepath, "centerline_on_pet.nii.gz"));
            cl_on_pet_flirt = copy(this.tof_on_pet_flirt);
            cl_on_pet_flirt.in = cl;
            cl_on_pet_flirt.out = this.centerline_on_pet_;
            cl_on_pet_flirt.interp = 'trilinear';
            cl_on_pet_flirt.applyXfm;
            assert(isfile(this.centerline_on_pet_.fqfn))   
            this.centerline_on_pet_.selectImagingTool();
            thr = dipmax(this.centerline_on_pet_)/10;
            this.centerline_on_pet_ = this.centerline_on_pet_.thresh(thr);
            this.centerline_on_pet_ = this.centerline_on_pet_.binarized();
            this.centerline_on_pet_.fileprefix = "centerline_on_pet";
            this.centerline_on_pet_.save();

            popd(pwd1);
        end
        function idif_ic = build_aif(this, petobj, tofobj, opts)
            arguments
                this mlaif.MipIdif
                petobj
                tofobj
                opts.cl = this.centerline_on_pet
                opts.mipt_thr double = 180000
                opts.frac_select double = 0.05
            end

            % e.g., !ls *tof*T1w*.nii.gz -> 'sub-108293_ses-20210218092914_tof_fl3d_tra_p2_multi-slab_orient-std_on_T1w.nii.gz'
            %tof__ = mlfourd.ImagingContext2(tofobj);

            % e.g., !ls *trc-oc*.nii.gz -> 'sub-108293_ses-20210421144815_trc-oc_proc-dyn_pet_on_T1w.nii.gz'
            tr = mlfourd.ImagingContext2(petobj);
            tr_single = single(tr);
            sz = size(tr_single);
            [tr_mipt,tr_indices] = max(tr, [], 4);
            cache.tr_mipt = tr_mipt;
            cache.tr_indices = tr_indices;
            %tr_mipt.view(opts.cl, tof__)

            % parse tracobj for tags
            re = mlpipeline.Bids.regexp_fileprefix(tr);
            assert(~isempty(re))
            if strcmp(re.trc, "co") || strcmp(re.trc, "oc")
                opts.mipt_thr = min(50000, opts.mipt_thr);
            end

            % select 0.05*numel of brightest and earliest-arriving voxels from centerline        
            centerline = logical(opts.cl);
            cache.centerline = centerline;

            tr_single__ = zeros(dipsum(centerline), sz(4));
            for t = 1:sz(4)
                vol_ = tr_single(:,:,:,t);
                tr_single__(:,t) = vol_(centerline);
            end
            cache.tr_single__ = tr_single__;

            brightest = single(tr_mipt);
            brightest = brightest(centerline);
            cache.brightest = brightest;

            earliest = single(tr_indices);
            earliest = earliest(centerline);
            cache.earliest = earliest;

            len = length(earliest);
            T = table(ascol(1:len), ascol(brightest), ascol(earliest), variableNames=["indices", "brightest", "earliest"]);
            T = sortrows(T, [2, 3], {'descend', 'ascend'});
            cache.T = T;
            T1 = T(ascol(1:round(opts.frac_select*len)), :);
            tr_single__ = tr_single__(T1.indices,:);
            tr_single__ = mean(tr_single__, 1);

            % save cache
            cache_fqfn = sprintf("%s_%s_cache.mat", this.pet_dyn.fqfp, stackstr());
            save(cache_fqfn, "cache");
            
            % aif from volume-average of brightest 0.05 from centerline
            idif_ic = this.save_imaging_context(mlfourd.ImagingContext2(tr_single__));

            % pcshow
            hold("on")
            toshow = tr_mipt.thresh(opts.mipt_thr); toshow.pcshow()
            opts.cl.pcshow()
            hold("off")
            saveFigure2(gcf, fullfile(tr.filepath, ...
                sprintf("pcshow_%s_%s_trc-%s_proc-tof-mips.fig", re.sub, re.ses, re.trc)))

            % plot
            h = figure;
            timesMid = idif_ic.json_metadata.timesMid; % previously made numeric
            plot(asrow(timesMid), asrow(idif_ic.imagingFormat.img), ':o')
            fontsize(scale=1.3)
            xlabel("times (s)")
            ylabel("IDIF activity density (Bq/mL)")
            title(sprintf("IDIF by MIPs, %s_%s_trc-%s", re.sub, re.ses, re.trc), Interpreter="none")
            savefig(h, strcat(idif_ic.fileprefix, ".fig"))
        end
        function idif_ic = build_all(this, opts)
            arguments
                this mlaif.MipIdif
                opts.steps logical = true(1, 5)
                opts.delete_large_files logical = true;
            end

            % cached on filesystem
            if isfile(this.new_fqfp + ".nii.gz")
                idif_ic = mlfourd.ImagingContext2(this.new_fqfp + ".nii.gz");
                return
            end

            if ~isfile(this.centerline_on_pet.fqfn)
                if opts.steps(1)
                    this.build_pet_objects(); % lots of fsl ops
                end
                if opts.steps(2)
                    this.build_tof_mips(); % manual drawing
                end
                if opts.steps(3)
                    this.back_project();
                end
            end
            if opts.steps(4)
                idif_ic = this.build_aif(this.pet_dyn, this.tof_on_pet); % by statistical sampling
            end
            if opts.steps(5) && ~strcmpi(this.tracer, "co") && ~strcmpi(this.tracer, "oc")
                idif_ic = this.build_deconv(idif_ic);
            end

            if opts.delete_large_files && ~contains(this.pet_dyn.filepath, "sourcedata")
                deleteExisting(this.pet_dyn.fqfp+".*")
            end

            %this.arterial_centerline_ = [];
            %this.arterial_input_function_ = [];
        end    
        function idif_ic = build_deconv(this, idif_ic)
            arguments 
                this mlaif.MipIdif
                idif_ic mlfourd.ImagingContext2
            end

            % deconvBayes
            boxcar = mlsiemens.BoxcarModel.create(artery=idif_ic, tracer=this.tracer); 
            boxcar.build_solution();
            idif_ic = boxcar.artery;
            idif_ic.fileprefix = mlpipeline.Bids.adjust_fileprefix(idif_ic.fileprefix, ...
                post_proc="deconv", remove_substring="_finite");
            idif_ic.save();
            %figure; plot(idif_ic)
        end
        function build_pet_objects(this)
            %% Build PET mips for use as guides when calling build_tof_mips.
            %  Also build transformations for registering centerlines to native PET.

            %% Register PET avgt on T1w.

            pet_on_t1w_fqfn = sprintf("%s_on_T1w.nii.gz", this.pet_avgt.fqfp);
            pet_on_t1w = mlfourd.ImagingContext2(pet_on_t1w_fqfn);
            pet_on_t1w_flirt = mlfsl.Flirt( ...
                'in', this.pet_avgt, ...
                'ref', this.t1w, ...
                'out', pet_on_t1w, ...
                'omat', this.mat(pet_on_t1w), ...
                'bins', 256, ...
                'cost', 'mutualinfo', ...
                'dof', 6, ...
                'interp', 'spline', ...
                'noclobber', false);
            if ~isfile(pet_on_t1w.fqfn)
                % do expensive coreg.
                pet_on_t1w_flirt.flirt();
            end
            assert(isfile(pet_on_t1w.fqfn))

            %% Register T1w on tof.

            t1w_on_tof = mlfourd.ImagingContext2(this.t1w.fqfp+"_on_tof.nii.gz");
            t1w_on_tof_flirt = mlfsl.Flirt( ...
                'in', this.t1w, ...
                'ref', this.tof, ...
                'out', t1w_on_tof, ...
                'omat', this.mat(t1w_on_tof), ...
                'bins', 256, ...
                'cost', 'mutualinfo', ...
                'dof', 6, ...
                'interp', 'spline', ...
                'noclobber', false);
            if ~isfile(t1w_on_tof.fqfn)
                % do expensive coreg.
                t1w_on_tof_flirt.flirt();
            end
            assert(isfile(t1w_on_tof.fqfn))

            %% Transform PET avgt to tof.

            this.pet_avgt_on_tof_ = mlfourd.ImagingContext2(this.pet_avgt.fqfp+"_on_tof.nii.gz");
            pet_avgt_on_tof_flirt = copy(pet_on_t1w_flirt);
            pet_avgt_on_tof_flirt.concatXfm(BtoC=this.mat(t1w_on_tof));
            pet_avgt_on_tof_flirt.in = this.pet_avgt;
            pet_avgt_on_tof_flirt.ref = this.tof;
            pet_avgt_on_tof_flirt.out = this.pet_avgt_on_tof_;
            if ~isfile(this.pet_avgt_on_tof_.fqfn)
                pet_avgt_on_tof_flirt.applyXfm();
            end
            assert(isfile(this.pet_avgt_on_tof_.fqfn))

            %% Transform PET mipt to tof.

            this.pet_mipt_on_tof_ = mlfourd.ImagingContext2(this.pet_mipt.fqfp+"_on_tof.nii.gz");
            pet_mipt_on_tof_flirt = copy(pet_avgt_on_tof_flirt);
            pet_mipt_on_tof_flirt.in = this.pet_mipt;
            pet_mipt_on_tof_flirt.out = this.pet_mipt_on_tof_;
            if ~isfile(this.pet_mipt_on_tof_.fqfn)
                pet_mipt_on_tof_flirt.applyXfm();
            end
            assert(isfile(this.pet_mipt_on_tof_.fqfn))

            %% Transform tof to PET mipt.

            fn = this.tof.fileprefix+"_on_pet_avgt.nii.gz";
            fqfn = fullfile(this.pet_dyn.filepath, fn);
            this.tof_on_pet_ = mlfourd.ImagingContext2(fqfn);
            this.tof_on_pet_flirt_ = copy(pet_avgt_on_tof_flirt);
            this.tof_on_pet_flirt_.invertXfm();
            this.tof_on_pet_flirt_.in = this.tof;
            this.tof_on_pet_flirt_.ref = this.pet_avgt;
            this.tof_on_pet_flirt_.out = this.tof_on_pet_;
            if ~isfile(this.tof_on_pet_) 
                this.tof_on_pet_flirt_.applyXfm();
            end
            assert(isfile(this.tof_on_pet_.fqfn))            
        end
        function build_tof_mips(this)
            
            pet_miptx = max(this.pet_mipt_on_tof, [], 1);
            pet_mipty = max(this.pet_mipt_on_tof, [], 2);
            pet_miptz = max(this.pet_mipt_on_tof, [], 3);

            %% draw in fsleyes 2-3 voxel widths for barycentric continuity

            pwd0 = pushd(this.pet_dyn.filepath);

            fprintf("Please draw arterial centerlines as overlaid images on the tof\n" + ...
                "Save files: centerline_zl and centerline_zr.\n\n");
            tof_mipz = max(this.tof, [], 3); 
            tof_mipz.view(pet_miptz) % pet_miptz.view(tof_mipz)
            fprintf("Please draw arterial centerlines as overlaid images on the tof\n" + ...
                "Save files: centerline_yl and centerline_yr.\n\n");
            tof_mipy = max(this.tof, [], 2); 
            tof_mipy.view(pet_mipty) % pet_mipty.view(tof_mipy)
            fprintf("Please draw arterial centerlines as overlaid images on the tof\n" + ...
                "Save files: centerline_xl and centerline_xr.\n\n");
            tof_mipx = max(this.tof, [], 1); 
            tof_mipx.view(pet_miptx) % pet_miptx.view(tof_mipx)

            popd(pwd0);
        end
        function out_ic = flirt_centerline(this, target, target_trc)
            arguments
                this mlkinetics.MipIdif
                target {mustBeNonempty}
                target_trc {mustBeTextScalar} = "unknown"
            end
            target = mlfourd.ImagingContext2(target);

            % flirt
            cl = this.centerline_on_pet;
            source = this.centerline_source;
            source_on_target = mlfourd.ImagingContext2( ...
                source.fqfp+"_on_"+target_trc+".nii.gz");
            cl_flirt = mlfsl.Flirt( ...
                'in', source, ...
                'ref', target, ...
                'out', source_on_target, ...
                'omat', this.mat(source_on_target), ...
                'bins', 256, ...
                'cost', 'mutualinfo', ...
                'dof', 6, ...
                'interp', 'spline', ...
                'noclobber', false);
            if ~isfile(source_on_trc.fqfn)
                % do expensive coreg.
                cl_flirt.flirt();
            end
            assert(isfile(source_on_trc.fqfn))

            % applyxfm
            out_ic = mlfourd.ImagingContext2( ...
                fullfile(target.filepath, "centerline_on_"+target_trc+"nii.gz"));
            cl_flirt.in = cl;
            cl_flirt.out = out_ic;
            cl_flirt.interp = "nearestneighbour";
            cl_flirt.applyXfm();
        end        
        function ff = new_fqfp(this, opts)
            %% seeks out pertinent fqfp in derivatives, not sourcedata

            arguments
                this mlaif.MipIdif
                opts.remove_substring {mustBeTextScalar} = "_timeAppend-4"
            end

            assert(contains(this.pet_dyn.filepath, "derivatives"))
            assert(~contains(this.pet_dyn.filepath, "sourcedata"))
            fp = mlpipeline.Bids.adjust_fileprefix(this.pet_dyn.fileprefix, ...                
                new_proc="MipIdif-finite-deconv", new_mode="idif", remove_substring=opts.remove_substring);
            ff = fullfile(this.pet_dyn.filepath, fp);
            if contains(ff, "*")
                ff = mglob(ff);
                assert(~isempty(ff))
                if numel(ff) > 1
                    warning("mlaif:RuntimeWarning", stackstr()+" returned string array of length "+numel(ff))
                end
            end
        end
        function g = new_fqfileprefix(this, varargin)
            g = this.new_fqfp(varargin{:});
        end
        function [ic] = pet_static_on_anatomy(this)
            ic = this.pet_avgt_on_tof;
        end
        function pet_dyn_on_anatomy(~)
            error("mlaif:NotImplementedError", "MipIdif");
        end
        function idif_ic = save_imaging_context(this, idif, opts)
            %% adjusts fileprefix and json_metadata prior to saving

            arguments
                this mlaif.MipIdif
                idif mlfourd.ImagingContext2
                opts.filepath = this.pet_dyn.filepath
                opts.fileprefix_template = this.pet_dyn.fileprefix
                opts.json_metadata_template = this.pet_dyn.json_metadata
            end

            idif_ic = mlfourd.ImagingContext2(idif);
            idif_ic.filepath = opts.filepath;

            % adjust fileprefix from template
            fp = opts.fileprefix_template;
            idif_ic.fileprefix = mlpipeline.Bids.adjust_fileprefix(fp, ...                
                new_proc="MipIdif", new_mode="idif", remove_substring="_timeAppend-4");

            % adjust json_metadata from template
            j = opts.json_metadata_template;
            idif_ic.json_metadata = j;
            idif_ic.selectNiftiTool();
            idif_ic.save();
            this.product_ = idif_ic;
        end
        function [aife,aifl] = splice_early_to_late(this, earlyobj, lateobj, len_late, trc, ident)
            arguments
                this mlaif.MipIdif
                earlyobj = "sub-108293_ses-20210421134537_fdg-dt0sec-BMC-LM-00-dynamic-3_interleave5_on_T1w.nii.gz"
                lateobj = "sub-108293_ses-20210421155709_trc-fdg_proc-dyn_pet_on_T1w.nii.gz"
                len_late = 19 % use last 19 frames only
                trc = "fdg"
                ident = "sub-108293_ses-20210421155709_trc-fdg"
            end
            
            M = 5; % # interleaves, which undercounts activity density
            
            cl = mlfourd.ImagingContext2('centerline_on_T1w.nii.gz');
            sz = size(cl);
            
            early = mlfourd.ImagingFormatContext2(earlyobj); % 0:2:308 sec
            sze = size(early);
            lene = sze(4) - 4;
            late = mlfourd.ImagingFormatContext2(lateobj); % 0:5:120, 120:20:300, 300:60:900, 900:300:3600
            lenl = len_late; 
            img = zeros(sz(1), sz(2), sz(3), lene+lenl);
            img(:,:,:,1:lene) = M*early.img(:,:,:,1:lene);
            img(:,:,:,lene+1:end) = late.img(:,:,:,end-lenl+1:end);
            early.img = img; % overwrite largest object
            clear img
            early.fileprefix = strrep(late.fileprefix, "proc-dyn", "proc-dyn-bmc-dt2tau10");
            early.save()
            
            early = mlfourd.ImagingContext2(early);
            aife = early.volumeAveraged(cl);
            aife.fileprefix = sprintf("aif_trc-%s_proc-tof-mips-early", trc);
            aife.save
            tause = [10*ones(1,146) 60*ones(1,10) 300*ones(1,9)]; % 1:146 interleaved with dt = 2 sec
            timesMide = nan(1,165);
            timesMide(1:146) = 5 + (0:2:290);
            timesMide(147:end) = 300 + cumsum(tause(147:end)) - tause(147:end)/2;
            
            plot(timesMide, aife.imagingFormat.img)
            xlabel("times (s)")
            ylabel("AIF activity (Bq/mL)")
            title(sprintf("AIF by TOF MIPs, %s", ident), Interpreter="none")
            
            late = mlfourd.ImagingContext2(late);
            aifl = late.volumeAveraged(cl);
            aifl.fileprefix = sprintf("aif_trc-%s_proc-tof-mips-late", trc);
            aifl.save
            tausl = [5*ones(1,24) 20*ones(1,9) 60*ones(1,10) 300*ones(1,9)];
            timesMidl = cumsum(tausl) - tausl/2;
            
            plot(timesMidl, aifl.imagingFormat.img)
            xlabel("times (s)")
            ylabel("AIF activity (Bq/mL)")
            title(sprintf("AIF by TOF MIPs, %s", ident), Interpreter="none")            
        end
    end

    methods (Static)
        function this = create(opts)
             arguments
                opts.bids_kit mlkinetics.BidsKit
                opts.tracer_kit mlkinetics.TracerKit
                opts.scanner_kit mlkinetics.ScannerKit
                opts.pet_avgt = []
                opts.pet_mipt = []
            end

            this = mlaif.MipIdif();
            this.bids_kit_ = opts.bids_kit;
            this.tracer_kit_ = opts.tracer_kit;
            this.scanner_kit_ = opts.scanner_kit;
            pet_dyn = this.scanner_kit_.do_make_activity_density(); % decayCorrected=true
            med = this.bids_kit_.make_bids_med();
            this.t1w_ = med.t1w_ic;
            this.tof_ = med.tof_ic;

            pet_dyn.relocateToDerivativesFolder();
            if ~isfile(pet_dyn)
                pet_dyn.ensureSingle;
                save(pet_dyn);
            end
            this.pet_dyn_fqfn_ = pet_dyn.fqfn;

            if ~isempty(opts.pet_avgt)
                this.pet_avgt_ = mlfourd.ImagingContext2(opts.pet_avgt);
            end
            if ~isempty(opts.pet_mipt)
                this.pet_mipt_ = mlfourd.ImagingContext2(opts.pet_mipt);
            end
        end   
        function fn = json(obj)
            if isa(obj, 'mlfourd.ImagingContext2')
                fn = strcat(obj.fqfp, '.json');
                return
            end
            ic = mlfourd.ImagingContext2(obj);
            fn = strcat(ic.fqfp, '.json');
        end
        function fn = mat(obj)
            if isa(obj, 'mlfourd.ImagingContext2')
                fn = strcat(obj.fqfp, '.mat');
                return
            end
            ic = mlfourd.ImagingContext2(obj);
            fn = strcat(ic.fqfp, '.mat');
        end
        function fn = niigz(obj)
            if isa(obj, 'mlfourd.ImagingContext2')
                fn = strcat(obj.fqfp, '.nii.gz');
                return
            end
            ic = mlfourd.ImagingContext2(obj);
            fn = strcat(ic.fqfp, '.nii.gz');
        end
        function viewAnat()
        end
    end

    %% PROTECTED

    properties (Access = protected)
        centerline_on_pet_
        pet_avgt_
        pet_avgt_on_tof_
        pet_dyn_fqfn_
        pet_mipt_
        pet_mipt_on_tof_
        product_
        t1w_
        tof_
        tof_on_pet_
        tof_on_pet_flirt_
        tracer_
        weights_timesMid_

        bids_kit_
        tracer_kit_
        scanner_kit_
    end

    methods (Access = protected)
        function this = MipIdif()
        end
        function that = copyElement(this)
            %%  See also web(fullfile(docroot, 'matlab/ref/matlab.mixin.copyable-class.html'))
            
            that = copyElement@matlab.mixin.Copyable(this);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
