classdef MipIdif < handle & mlsystem.IHandle
    %% This is a builder for generating centerlines, then applying them to dynamic PET to obtain IDIFs.  
    %  Start with build_all(this).  
    %  
    %  Created 12-Jun-2023 22:59:36 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlaif/src/+mlaif.
    %  Developed on Matlab 9.14.0.2254940 (R2023a) Update 2 for MACI64.  Copyright 2023 John J. Lee.
    
    properties (Constant)
        ALPHA = 0.25
        FLIRT_COST = 'normmi'
        max_len_mipt = 30
    end

    properties
        filename_pattern_dynamic
        filename_pattern_dynamic_fdg
        filename_pattern_static
        filename_pattern_static_fdg
        minz_for_mip
        reference_tracer
    end
    
    properties (Dependent)
        centerline_on_reftrc
        centerline_on_pet
        fdg_dyn
        fdg_avgt
        filename_pattern_dynamic_ref
        filename_pattern_static_ref
        halflife
        json_metadata_template
        model_kind
        pet_avgt % w/ decay-correction, weighted by taus
        pet_dyn % read from filesystem to minimize copies of objects in memory, relocated from sourcedata -> derivatives
        pet_mipt % w/ decay-correction
        product
        reftrc_dyn
        reftrc_avgt
        t1w
        t1w_on_pet
        t1w_on_pet_fqfn % not necessarily a file
        tof
        tof_on_pet
        tof_on_pet_fqfn % not necessarily a file
        tof_on_t1w
        tracer
        weights_timesMid
    end

    methods %% GET, SET
        function g = get.centerline_on_reftrc(this)
            if isempty(this.centerline_on_reftrc_)
                try
                    pth = this.reftrc_avgt.filepath;
                    this.centerline_on_reftrc_ = mlfourd.ImagingContext2( ...
                        fullfile(pth, "centerline_on_pet_MipIdif.nii.gz"));
                catch ME
                    handexcept(ME)
                end
            end
            g = copy(this.centerline_on_reftrc_);
        end
        function g = get.centerline_on_pet(this)
            if isempty(this.centerline_on_pet_)
                try
                    if this.dilate_m_ > 0
                        this.centerline_on_pet_ = mlfourd.ImagingContext2( ...
                            fullfile(this.pet_mipt.filepath, ...
                            sprintf("centerline_on_pet_dilate-m%i.nii.gz", this.dilate_m_)));
                    else
                        this.centerline_on_pet_ = mlfourd.ImagingContext2( ...
                            fullfile(this.pet_mipt.filepath, "centerline_on_pet_MipIdif.nii.gz"));
                    end
                    % assert(isfile(this.centerline_on_pet_))
                catch ME
                    handexcept(ME)
                end
            end
            g = copy(this.centerline_on_pet_);
        end
        function g = get.fdg_dyn(this)
            if ~isempty(this.fdg_dyn_)
                g = this.fdg_dyn_;
                return
            end

            % dynamic
            mg = mglob(fullfile(this.mediator_.sourceSubPath, "ses-*", "pet", this.filename_pattern_dynamic_fdg));
            mg = natsort(mg); % human fdg always precedes phantom fdg; prefer delay0
            if ~isemptytext(mg) && isfile(mg(1))
                % pick delay0
                this.fdg_dyn_ = mlfourd.ImagingContext2(mg(1));
                % avoid relocating large dynamic imaging to derivatives
                % this.reftrc_avgt_.relocateToDerivativesFolder(); 
                g = this.fdg_dyn_;
                return
            end

            error("mlaif:RunTimeError", stackstr())
        end
        function g = get.fdg_avgt(this)
            if ~isempty(this.fdg_avgt_)
                g = this.fdg_avgt_;
                return
            end

            % prefer static
            mg = mglob(fullfile(this.mediator_.sourceSubPath, "ses-*", "pet", this.filename_pattern_static_fdg));
            mg = natsort(mg); % human fdg always precedes phantom fdg; prefer delay0
            if ~isemptytext(mg) && isfile(mg(1))
                % pick delay0
                this.fdg_avgt_ = mlfourd.ImagingContext2(mg(1));
                this.fdg_avgt_.relocateToDerivativesFolder();
                g = this.fdg_avgt_;
                return
            end

            error("mlaif:RunTimeError", stackstr())
        end
        function g = get.filename_pattern_dynamic_ref(this)
            g = "sub-*_ses-*_trc-"+this.reference_tracer+"_proc-delay*BrainMoCo2-createNiftiMovingAvgFrames.nii.gz";
        end
        function g = get.filename_pattern_static_ref(this)
            g = "sub-*_ses-*_trc-"+this.reference_tracer+"_proc-delay*BrainMoCo2-createNiftiStatic.nii.gz";
        end
        function g = get.halflife(this)
            rn = this.tracer_kit_.make_radionuclides();
            g = rn.halflife;
            assert(isfinite(g))
        end
        function g = get.json_metadata_template(this)
            if ~isempty(this.json_metadata_template_)
                g = this.json_metadata_template_;
                return
            end

            j = this.pet_dyn.json_metadata;
            this.json_metadata_template_ = mlpipeline.ImagingMediator.ensureNumericTimingData(j);
            g = this.json_metadata_template_;
        end
        function g = get.model_kind(this)
            g = this.model_kind_;
        end
        function g = get.pet_avgt(this)
            if ~isempty(this.pet_avgt_)
                g = this.pet_avgt_;
                return
            end

            if contains(this.tracer, "oo", IgnoreCase=true)
                % OO static is often corrupted by motion
                this.pet_avgt_ = this.pet_dyn.timeAveraged(weights=asrow(this.pet_dyn.json_metadata.taus));
                this.pet_avgt_.relocateToDerivativesFolder();
                if ~isfile(this.pet_avgt_)
                    save(this.pet_avgt_);
                end
                g = this.pet_avgt_;
            end

            % prefer static
            pth = fullfile(this.mediator_.sourceSesPath, "pet");
            mg = mglob(fullfile(pth, this.filename_pattern_static));
            if ~isemptytext(mg)
                % pick delay0
                select = contains(mg, "delay0");
                this.pet_avgt_ = mlfourd.ImagingContext2(mg(select));
                this.pet_avgt_.relocateToDerivativesFolder();
                g = this.pet_avgt_;
                return
            end

            % console recons
            mg = mglob(fullfile(pth, sprintf("sub-*_ses-*_trc-%s_proc-consoleStatic.nii.gz", this.tracer)));
            if ~isemptytext(mg)
                this.pet_avgt_ = mlfourd.ImagingContext2(mg(1));
                this.pet_avgt_.relocateToDerivativesFolder();
                g = this.pet_avgt_;
                return
            end

            % look for _avgt
            fqfn = this.pet_dyn.fqfp+"_avgt.nii.gz";
            if isfile(fqfn)
                this.pet_avgt_ = mlfourd.ImagingContext2(fqfn);
                this.pet_avgt_.relocateToDerivativesFolder();
                g = this.pet_avgt_;
                return
            end

            % construct _avgt
            this.pet_avgt_ = this.pet_dyn.timeAveraged(weights=asrow(this.pet_dyn.json_metadata.taus));
            this.pet_avgt_.relocateToDerivativesFolder();
            if ~isfile(this.pet_avgt_)
                save(this.pet_avgt_);
            end
            g = this.pet_avgt_;
        end
        function g = get.pet_dyn(this)
            if ~isempty(this.pet_dyn_)
                g = this.pet_dyn_;
                return
            end

            ic = this.mediator_.imagingContext;
            if contains(ic.fqfn, "-delay0") && ~contains(ic.fqfn, "_timeAppend")
                g = mlkinetics.BidsKit.timeAppend(ic, do_save=true);  % also manages trivial cases with nothing to append
            else 
                g = ic;
            end
            this.pet_dyn_ = g;
        end
        function g = get.pet_mipt(this)
            if ~isempty(this.pet_mipt_)
                g = this.pet_mipt_;
                return
            end

            fqfn = this.pet_dyn.fqfp+"_mipt.nii.gz";
            fqfn = strrep(fqfn, "sourcedata", "derivatives");
            if isfile(fqfn)
                this.pet_mipt_ = mlfourd.ImagingContext2(fqfn);
                g = this.pet_mipt_;
                return
            end

            ifc = this.pet_dyn.imagingFormat;   
            L = min(size(ifc, 4), this.max_len_mipt);
            ifc.img = ifc.img(:,:,:,1:L);
            ifc.img(:,:,1:this.minz_for_mip-1,:) = 0;
            ifc.img(:,:,end-this.minz_for_mip+1:end,:) = 0;
            % weights = this.weights_timesMid(1:L); % negligible differences
            % for t = 1:L
            %     ifc.img(:,:,:,t) = ifc.img(:,:,:,t)*weights(t);
            % end

            ic = mlfourd.ImagingContext2(ifc);
            this.pet_mipt_ = max(ic, [], 4);
            this.pet_mipt_.fileprefix = sprintf("%s_mipt", this.pet_dyn.fileprefix);
            this.pet_mipt_.relocateToDerivativesFolder();
            if ~isfile(this.pet_mipt_)
                save(this.pet_mipt_);
            end
            g = this.pet_mipt_;
        end
        function g = get.product(this)
            assert(~isempty(this.product_))
            g = copy(this.product_);
        end
        function g = get.reftrc_dyn(this)
            if ~isempty(this.reftrc_dyn_)
                g = this.reftrc_dyn_;
                return
            end

            % dynamic
            mg = mglob(fullfile(this.mediator_.sourceSubPath, "ses-*", "pet", this.filename_pattern_dynamic_ref));
            mg = natsort(mg); % human fdg always precedes phantom fdg; prefer delay0
            if ~isemptytext(mg) && isfile(mg(1))
                % pick delay0
                this.reftrc_dyn_ = mlfourd.ImagingContext2(mg(1));
                % avoid relocating large dynamic imaging to derivatives
                % this.reftrc_avgt_.relocateToDerivativesFolder(); 
                g = this.reftrc_dyn_;
                return
            end

            error("mlaif:RunTimeError", stackstr())
        end
        function g = get.reftrc_avgt(this)
            if ~isempty(this.reftrc_avgt_)
                g = this.reftrc_avgt_;
                return
            end

            % prefer static
            mg = mglob(fullfile(this.mediator_.sourceSubPath, "ses-*", "pet", this.filename_pattern_static_ref));
            mg = natsort(mg); % human fdg always precedes phantom fdg; prefer delay0
            if ~isemptytext(mg) && isfile(mg(1))
                % pick delay0
                this.reftrc_avgt_ = mlfourd.ImagingContext2(mg(1));
                this.reftrc_avgt_.relocateToDerivativesFolder();
                g = this.reftrc_avgt_;
                return
            end

            error("mlaif:RunTimeError", stackstr())
        end
        function g = get.t1w(this)
            g = this.mediator_.t1w_ic;
        end
        function g = get.t1w_on_pet(this)
            if ~isempty(this.t1w_on_pet_)
                g = this.t1w_on_pet_;
                return
            end

            if ~isempty(this.t1w_on_pet_fqfn) && isfile(this.t1w_on_pet_fqfn)
                this.t1w_on_pet_ = mlfourd.ImagingContext2(this.t1w_on_pet_fqfn);
                g = this.t1w_on_pet_;
                return
            end

            g = [];
        end
        function g = get.t1w_on_pet_fqfn(this)
            %% not necessarily a file

            if ~isempty(this.t1w_on_pet_fqfn_)
                g = this.t1w_on_pet_fqfn_;
                return
            end

            this.t1w_on_pet_fqfn_ = fullfile(this.pet_avgt.filepath, "T1w_on_"+this.pet_avgt.fileprefix+".nii.gz");
            this.t1w_on_pet_fqfn_ = strrep(this.t1w_on_pet_fqfn_, "sourcedata", "derivatives");
            g = this.t1w_on_pet_fqfn_;
        end
        function g = get.tof(this)
            g = this.mediator_.tof_ic;
        end
        function g = get.tof_on_pet(this)
            if ~isempty(this.tof_on_pet_)
                g = this.tof_on_pet_;
                return
            end

            if ~isempty(this.tof_on_pet_fqfn) && isfile(this.tof_on_pet_fqfn)
                this.tof_on_pet_ = mlfourd.ImagingContext2(this.tof_on_pet_fqfn);
                g = this.tof_on_pet_;
                return
            end

            g = [];
        end
        function g = get.tof_on_pet_fqfn(this)
            %% not necessarily a file

            if ~isempty(this.tof_on_pet_fqfn_)
                g = this.tof_on_pet_fqfn_;
                return
            end

            this.tof_on_pet_fqfn_ = fullfile(this.pet_avgt.filepath, "tof_on_"+this.pet_avgt.fileprefix+".nii.gz");
            this.tof_on_pet_fqfn_ = strrep(this.tof_on_pet_fqfn_, "sourcedata", "derivatives");
            g = this.tof_on_pet_fqfn_;
        end
        function g = get.tof_on_t1w(this)
            g = this.mediator_.tof_on_t1w_ic;
        end
        function g = get.tracer(this)
            if ~isempty(this.tracer_)
                g = this.tracer_;
                return
            end

            g = this.mediator_.tracer;
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
        function cl = align_centerline_on_tracer(this, opts)
            %% aligns centerline from reference tracer on other tracers

            arguments
                this mlaif.MipIdif
                opts.frames double = []
                opts.stat function_handle = @max
                opts.do_apply_mask logical = false
                opts.blur_for_mask {mustBeScalarOrEmpty} = 20  % mm blurring of mask at fwhh
                opts.interp {mustBeTextScalar} = "nearestneighbour"  % flirt | grep interpolation
                opts.noclobber logical = false
            end

            if contains(this.pet_dyn.fqfn, "_trc-"+this.reference_tracer, IgnoreCase=true)  % ref needs no aligning
                cl = this.centerline_on_pet;
                assert(isfile(cl));
                return
            end

            % use early ref frames as needed to align vascular geometry
            if ~isempty(opts.frames)
                reftrac_avgt__ = this.select_frames(this.reftrc_dyn, frames=opts.frames, stat=opts.stat);
            else
                reftrac_avgt__ = this.reftrc_avgt;
            end

            % use dlicv masks
            inweight = this.build_dlicv_mask(this.pet_avgt, blur=opts.blur_for_mask);
            refweight = this.build_dlicv_mask(reftrac_avgt__, blur=opts.blur_for_mask);

            % flirt this.pet_avgt to reftrac_avgt__
            pet_on_ref_fqfn = this.pet_avgt.fqfp+"_on_"+this.reference_tracer+".nii.gz"; 
            pet_on_ref = mlfourd.ImagingContext2(pet_on_ref_fqfn);
            if opts.do_apply_mask            
                pet_on_ref_flirt = mlfsl.Flirt( ...
                    'in', this.apply_mask(this.pet_avgt, inweight), ...
                    'inweight', inweight, ...
                    'ref', this.apply_mask(reftrac_avgt__, refweight), ...
                    'refweight', refweight, ...
                    'out', pet_on_ref, ...
                    'omat', this.mat(pet_on_ref), ...
                    'bins', 256, ...
                    'cost', this.FLIRT_COST, ...
                    'dof', 6, ...
                    'interp', 'spline', ...
                    'noclobber', false);
            else
                pet_on_ref_flirt = mlfsl.Flirt( ...
                    'in', this.pet_avgt, ...
                    'inweight', inweight, ...
                    'ref', reftrac_avgt__, ...
                    'refweight', refweight, ...
                    'out', pet_on_ref, ...
                    'omat', this.mat(pet_on_ref), ...
                    'bins', 256, ...
                    'cost', this.FLIRT_COST, ...
                    'dof', 6, ...
                    'interp', 'spline', ...
                    'noclobber', false);
            end
            if ~opts.noclobber || ~isfile(pet_on_ref.fqfn)
                % do expensive coreg.
                pet_on_ref_flirt.flirt();
            end
            assert(isfile(pet_on_ref.fqfn))

            % flirt this.centerline_on_reftrc to this.pet_avgt
            pet_on_ref_flirt.invertXfm();
            pet_on_ref_flirt.in = this.centerline_on_reftrc;
            pet_on_ref_flirt.out = this.centerline_on_pet;
            pet_on_ref_flirt.ref = this.pet_avgt;
            pet_on_ref_flirt.interp = opts.interp;
            pet_on_ref_flirt.applyXfm();
            copyfile(this.centerline_on_pet.fqfn, strrep(this.centerline_on_pet.fqfn, "MipIdif", opts.interp))
            if ~strcmp(opts.interp, "nearestneighbour")
                this.centerline_on_pet_ = this.centerline_on_pet.thresh(this.ALPHA);
                this.centerline_on_pet_ = this.centerline_on_pet.binarized();
            end
            this.centerline_on_pet_.fileprefix = "centerline_on_pet_MipIdif";
            this.centerline_on_pet_.save();
            assert(isfile(this.centerline_on_pet))
        end

        function cl = align_centerline_on_tracers(this, opts)
            %% aligns centerline from reference tracer on other tracers, using FDG intermediary;
            %  default parameters support HO centerline -> FDG -> CO

            arguments
                this mlaif.MipIdif
                opts.frames double = 1:5
                opts.stat function_handle = @max
                opts.do_apply_mask logical = false
                opts.blur_for_mask {mustBeScalarOrEmpty} = 20  % mm blurring of mask at fwhh
                opts.interp {mustBeTextScalar} = "nearestneighbour"  % flirt | grep interpolation
                opts.noclobber logical = false
            end

            if contains(this.pet_dyn.fqfn, "_trc-"+this.reference_tracer, IgnoreCase=true)  % ref needs no aligning
                cl = this.centerline_on_pet;
                assert(isfile(cl));
                return
            end

            % use early FDG frames as needed to align vascular geometry
            if ~isempty(opts.frames)
                fdg_early__ = this.select_frames(this.fdg_dyn, frames=opts.frames, stat=opts.stat);
            else
                fdg_early__ = this.fdg_avgt;
            end

            % use dlicv masks
            inweight = [];  % this.build_dlicv_mask(this.pet_avgt, blur=opts.blur_for_mask);
            refweight = [];  % this.build_dlicv_mask(fdg_early__, blur=opts.blur_for_mask);

            % flirt this.pet_avgt to fdg_early__
            pet_on_fdg_fqfn = this.pet_avgt.fqfp+"_on_fdg.nii.gz"; 
            pet_on_fdg = mlfourd.ImagingContext2(pet_on_fdg_fqfn);
            if opts.do_apply_mask            
                pet_on_fdg_flirt = mlfsl.Flirt( ...
                    'in', this.apply_mask(this.pet_avgt, inweight), ...
                    'inweight', inweight, ...
                    'ref', this.apply_mask(fdg_early__, refweight), ...
                    'refweight', refweight, ...
                    'out', pet_on_fdg, ...
                    'omat', this.mat(pet_on_fdg), ...
                    'bins', 256, ...
                    'cost', this.FLIRT_COST, ...
                    'dof', 6, ...
                    'interp', 'spline', ...
                    'noclobber', false);
            else
                pet_on_fdg_flirt = mlfsl.Flirt( ...
                    'in', this.pet_avgt, ...
                    'inweight', inweight, ...
                    'ref', fdg_early__, ...
                    'refweight', refweight, ...
                    'out', pet_on_fdg, ...
                    'omat', this.mat(pet_on_fdg), ...
                    'bins', 256, ...
                    'cost', this.FLIRT_COST, ...
                    'dof', 6, ...
                    'interp', 'spline', ...
                    'noclobber', false);
            end
            if ~opts.noclobber || ~isfile(pet_on_fdg.fqfn)
                % do expensive coreg.
                pet_on_fdg_flirt.flirt();
            end
            assert(isfile(pet_on_fdg.fqfn))

            % use dlicv masks
            inweight1 = this.build_dlicv_mask(this.fdg_avgt, blur=opts.blur_for_mask);
            refweight1 = this.build_dlicv_mask(this.reftrc_avgt, blur=opts.blur_for_mask);

            % flirt this.fdg_avgt to this.reftrc_avgt
            fdg_on_reftrc_fqfn = this.fdg_avgt.fqfp+"_on_"+this.reference_tracer+".nii.gz"; 
            fdg_on_reftrc = mlfourd.ImagingContext2(fdg_on_reftrc_fqfn);
            if opts.do_apply_mask            
                fdg_on_reftrc_flirt = mlfsl.Flirt( ...
                    'in', this.apply_mask(this.fdg_avgt, inweight1), ...
                    'inweight', inweight1, ...
                    'ref', this.apply_mask(this.reftrc_avgt, refweight1), ...
                    'refweight', refweight1, ...
                    'out', fdg_on_reftrc, ...
                    'omat', this.mat(fdg_on_reftrc), ...
                    'bins', 256, ...
                    'cost', this.FLIRT_COST, ...
                    'dof', 6, ...
                    'interp', 'spline', ...
                    'noclobber', false);
            else
                fdg_on_reftrc_flirt = mlfsl.Flirt( ...
                    'in', this.fdg_avgt, ...
                    'inweight', inweight1, ...
                    'ref', this.reftrc_avgt, ...
                    'refweight', refweight1, ...
                    'out', fdg_on_reftrc, ...
                    'omat', this.mat(fdg_on_reftrc), ...
                    'bins', 256, ...
                    'cost', this.FLIRT_COST, ...
                    'dof', 6, ...
                    'interp', 'spline', ...
                    'noclobber', false);
            end
            if ~opts.noclobber || ~isfile(fdg_on_reftrc.fqfn)
                % do expensive coreg.
                fdg_on_reftrc_flirt.flirt();
            end
            assert(isfile(fdg_on_reftrc.fqfn))

            % concat transformations:  this.pet_avgt -> this.fdg_avgt > this.reftrc_avgt
            pet_avgt_on_reftrc_flirt = copy(fdg_on_reftrc_flirt);
            pet_avgt_on_reftrc_flirt.concatXfm(AtoB=this.mat(pet_on_fdg));

            % invert transformations:  this.reftrc_avgt -> this.fdg_avgt -> this.pet_avgt
            pet_avgt_on_reftrc_flirt.invertXfm();
            pet_avgt_on_reftrc_flirt.in = this.centerline_on_reftrc;
            pet_avgt_on_reftrc_flirt.out = this.centerline_on_pet;
            pet_avgt_on_reftrc_flirt.ref = this.pet_avgt;
            pet_avgt_on_reftrc_flirt.interp = opts.interp;
            pet_avgt_on_reftrc_flirt.applyXfm();
            copyfile(this.centerline_on_pet.fqfn, strrep(this.centerline_on_pet.fqfn, "MipIdif", opts.interp))
            if ~strcmp(opts.interp, "nearestneighbour")
                this.centerline_on_pet_ = this.centerline_on_pet.thresh(this.ALPHA);
                this.centerline_on_pet_ = this.centerline_on_pet.binarized();
            end
            this.centerline_on_pet_.fileprefix = "centerline_on_pet_MipIdif";
            this.centerline_on_pet_.save();
            assert(isfile(this.centerline_on_pet))
        end
        
        function ic = build_dlicv_mask(this, pet_static, opts)
            arguments
                this mlaif.MipIdif
                pet_static {mustBeNonempty}
                opts.blur {mustBeScalarOrEmpty} = []  % 5 mm fwhh blur; flirt fails with 40 mm
            end
            pet_static = mlfourd.ImagingContext2(pet_static);

            t1w_on_pet_flirt = mlfsl.Flirt( ...
                'in', this.t1w, ...
                'ref', pet_static, ...
                'out', this.t1w_on_pet, ...
                'omat', this.mat(this.t1w_on_pet), ...
                'bins', 256, ...
                'cost', 'mutualinfo', ...
                'dof', 6, ...
                'interp', 'spline', ...
                'noclobber', false);
            if ~isfile(this.t1w_on_pet_fqfn) || ~isfile(this.mat(this.t1w_on_pet))
                % do expensive coreg.
                t1w_on_pet_flirt.flirt();
            end

            ic_fqfn = pet_static.fqfp+"_dlicv.nii.gz";
            ic = mlfourd.ImagingContext2(ic_fqfn);
            if ~isfile(ic.fqfn)
                t1w_on_pet_flirt.in = this.mediator_.dlicv_ic;
                t1w_on_pet_flirt.out = ic;
                t1w_on_pet_flirt.interp = "nearestneighbour";
                t1w_on_pet_flirt.applyXfm();
            end
            assert(isfile(ic.fqfn))

            if ~isempty(opts.blur) && opts.blur > 0
                ic = ic.blurred(opts.blur);
            end
        end
        
        function cl = back_project(this)
            %% back-project MIPs            
            
            pwd1 = pushd(this.pet_avgt.filepath);

            % read_avw <- $FSLDIR/etc/matlab
            cxL_img = read_avw(fullfile(this.pet_avgt.filepath, "centerline_xl.nii.gz"));
            cxR_img = read_avw(fullfile(this.pet_avgt.filepath, "centerline_xr.nii.gz"));
            cyL_img = read_avw(fullfile(this.pet_avgt.filepath, "centerline_yl.nii.gz"));
            cyR_img = read_avw(fullfile(this.pet_avgt.filepath, "centerline_yr.nii.gz"));
            czL_img = read_avw(fullfile(this.pet_avgt.filepath, "centerline_zl.nii.gz"));
            czR_img = read_avw(fullfile(this.pet_avgt.filepath, "centerline_zr.nii.gz"));
            
            size_pet_avgt = size(this.pet_avgt);
            img1 = zeros(size_pet_avgt);
            img2 = zeros(size_pet_avgt);
            img3 = zeros(size_pet_avgt);
            img4 = zeros(size_pet_avgt);
            img5 = zeros(size_pet_avgt);
            img6 = zeros(size_pet_avgt);
            
            for ix = 1:size_pet_avgt(1)
                img1(ix,:,:) = cxL_img; end
            for iy = 1:size_pet_avgt(2)
                img2(:,iy,:) = cyL_img; end
            for iz = 1:size_pet_avgt(3)
                img3(:,:,iz) = czL_img; end
            for ix = 1:size_pet_avgt(1)
                img4(ix,:,:) = cxR_img; end
            for iy = 1:size_pet_avgt(2)
                img5(:,iy,:) = cyR_img; end
            for iz = 1:size_pet_avgt(3)
                img6(:,:,iz) = czR_img; end
            
            % img1 = imdilate(img1, strel("sphere", 1));
            % img2 = imdilate(img2, strel("sphere", 1));
            % img3 = imdilate(img3, strel("sphere", 1));
            % img4 = imdilate(img4, strel("sphere", 1));
            % img5 = imdilate(img5, strel("sphere", 1));
            % img6 = imdilate(img6, strel("sphere", 1));
            
            imgL = img1 & img2 & img3;
            imgR = img4 & img5 & img6;
            
            cl = copy(this.pet_avgt); 
            cl = cl.nifti;
            cl.img = imgL + imgR;
            cl = mlfourd.ImagingContext2(cl);
            cl = cl.binarized();
            %cl.pcshow()  
            cl.filepath = this.new_filepath;
            cl.fileprefix = 'centerline_on_pet_MipIdif';
            cl.save()
            this.centerline_on_pet_ = cl;

            popd(pwd1);
        end
        
        function idif_ic = build_aif(this, petobj, opts)
            arguments
                this mlaif.MipIdif
                petobj
                opts.cl = this.centerline_on_pet
                opts.mipt_thr double = 180000
                opts.frac_select double = this.ALPHA
                opts.do_view logical = false
            end

            cl_tags = split(opts.cl.fileprefix, "centerline_on_pet");
            cl_tags = strrep(cl_tags(2), "_", "-");
            sel_tags = sprintf("-select%g", opts.frac_select);
            sel_tags = strrep(sel_tags, ".", "p");
            tags = cl_tags + sel_tags;
            
            % adjust opts.mipt_thr for CO
            tr = mlfourd.ImagingContext2(petobj);
            re = mlpipeline.Bids.regexp_fileprefix(tr);
            assert(~isempty(re))
            if strcmp(re.trc, "co") || strcmp(re.trc, "oc")
                opts.mipt_thr = min(50000, opts.mipt_thr);
            end

            %% select this.ALPHA*numel of brightest and earliest-arriving voxels from centerline        

            centerline = logical(opts.cl);
            cache.centerline = centerline;

            % tr_mipt, tr_indices ~ MIP, time of max
            tr_single = single(tr);
            sz = size(tr_single);
            [tr_mipt,tr_indices] = max(tr, [], 4); % ~ 45 sec for dyn CO with 200 frames
            if ~isfile(tr_mipt)                    
                tr_mipt.save();
            end
            if ~isfile(tr_indices)
                tr_indices.save();
            end
            cache.tr_mipt = tr_mipt;
            cache.tr_indices = tr_indices;
            if opts.do_view
                tr_mipt.view();
                tr_indices.view();
            end

            % vec containing values of brightest 
            brightest = single(tr_mipt);
            brightest = brightest(centerline);
            cache.brightest = brightest;

            % vec containing index of earliest peaks
            earliest = single(tr_indices);
            earliest = earliest(centerline);
            cache.earliest = earliest;

            % tr_single__ ~ time-series selected from centerline
            tr_single__ = zeros(dipsum(centerline), sz(4));
            for t = 1:sz(4)
                vol_ = tr_single(:,:,:,t);
                tr_single__(:,t) = vol_(centerline);
            end
            cache.tr_single__ = tr_single__;

            % sort brightest, earliest
            len = length(earliest);
            T = table(ascol(1:len), ascol(brightest), ascol(earliest), variableNames=["indices", "brightest", "earliest"]);
            T = sortrows(T, [2, 3], {'descend', 'ascend'});
            cache.T = T;
            len_selected = max(1, round(opts.frac_select*len));
            T1 = T(ascol(1:len_selected), :);
            cache.T1 = T1;
            tr_single___ = tr_single__(T1.indices,:);
            cache.tr_single___ = tr_single___;

            % save cache
            cache_fqfn = sprintf("%s_%s%s_cache.mat", this.new_fqfp, stackstr(), tags);
            save(cache_fqfn, "cache");
            
            % aif from volume-average of brightest fraction from centerline ROI
            idif_ic = this.save_imaging_context( ...
                mlfourd.ImagingContext2(mean(tr_single___, 1)), ...
                tags=tags);

            if opts.do_view
                h = figure;
                timesMid = this.json_metadata_template.timesMid;
                for pidx = 1:size(T1, 1)
                    plot(asrow(timesMid), asrow(tr_single___(pidx, :)));
                end
                atitle = stackstr(use_dashes=true)+"-brightest-earliest";
                title(atitle)
                xlabel("times (s)")
                ylabel("activity (Bq/mL)")
                saveFigure2(h, ...
                    fullfile(opts.cl.filepath, atitle), ...
                    closeFigure=false);
            end

            if opts.do_view
                % pcshow
                h = figure;
                hold("on")
                toshow = tr_mipt.thresh(opts.mipt_thr); toshow.pcshow()
                opts.cl.pcshow()
                hold("off")
                saveFigure2(h, fullfile(tr.filepath, ...
                    sprintf("pcshow_%s_%s_trc-%s_proc-%s.fig", re.sub, re.ses, re.trc, stackstr(use_dashes=true))))

                % plot
                h1 = figure;
                timesMid = idif_ic.json_metadata.timesMid; % previously made numeric
                plot(asrow(timesMid), asrow(idif_ic.imagingFormat.img), ':o')
                %fontsize(scale=1.3)
                xlabel("times (s)")
                ylabel("IDIF activity density (Bq/mL)")
                title(sprintf("IDIF by MIPs, %s_%s_trc-%s", re.sub, re.ses, re.trc), Interpreter="none")
                savefig(h1, strcat(idif_ic.fileprefix, ".fig"))
            end
        end
        
        function idif_ic = build_all(this, opts)
            %% Args:
            %  opts.delete_large_files logical = false
            %  opts.steps logical = true(1, 6)
            %
            %       step(1) ~ calls build_pet_objects()
            %       step(2) ~ calls draw_tof_mips()
            %       step(3) ~ calls back_project()
            %       step(4) ~ calls align_centerline_on_tracer()
            %       step(5) ~ calls build_aif(this.pet_dyn)
            %       step(6) ~ calls build_deconv()

            arguments
                this mlaif.MipIdif
                opts.steps logical = true(1, 6)
                opts.delete_large_files logical = false
                opts.frac_select double = this.ALPHA
            end

            % cached on filesystem
            if isfile(this.new_fqfp + ".nii.gz")
                idif_ic = mlfourd.ImagingContext2(this.new_fqfp + ".nii.gz");
                return
            end

            idif_ic = [];
            if opts.steps(1)
                this.build_pet_objects(); % lots of fsl ops
            end
            if opts.steps(2)
                % try
                %     this.draw_tof_centerline();
                % catch ME
                %     fprintf("%s: %s\n", stackstr(), ME.message);  
                %     fprintf("%s: tof not available, trying T1w\n", stackstr());
                %     this.draw_t1w_centerline();
                % end
            end
            if opts.steps(3)
                % this.back_project(); % overwrites centerline_on_pet.nii.gz
            end
            if opts.steps(4)
                if contains(this.tracer, "co", IgnoreCase=true) || ...
                        contains(this.tracer, "oc", IgnoreCase=true)
                    if contains(this.reference_tracer, "fdg", IgnoreCase=true)
                        this.align_centerline_on_tracer(frames=1:25, stat=@max, do_apply_mask=false);
                    elseif contains(this.reference_tracer, "ho", IgnoreCase=true)
                        this.align_centerline_on_tracers(frames=1:6);
                    else
                        error("mlaif:RuntimeError", stackstr())
                    end
                else
                    this.align_centerline_on_tracer();
                end
            end
            if opts.steps(5)
                idif_ic = this.build_aif(this.pet_dyn, frac_select=opts.frac_select); % by statistical sampling
            end
            if opts.steps(6) % && ~strcmpi(this.tracer, "co") && ~strcmpi(this.tracer, "oc")
                % idif_ic = this.build_deconv(idif_ic);
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

            pwd0 = pushd(idif_ic.filepath);

            boxcar = mlnest.Boxcar.create(artery=idif_ic, model_kind=this.model_kind);
            obj = mlnest.MultiNest(context=boxcar);
            obj.filepath = idif_ic.filepath;
            obj = obj.solve( ...
                signal_model=@boxcar.signalmodel, ...
                verbose=false, ...
                Nlive=55, ...
                Nmcmc=0); 

            disp(obj)
            fprintf("Multinest.product: \n"); 
            disp(obj.product)
            fprintf("ks(): %g\n", obj.ks())
            fprintf("loss: %s\n", obj.loss())
            obj.plot_posteriors(singles=true, do_save=true);
            obj.saveall();

            [sig,idl] = obj.simulate();

            figure;
            plot_over_figure(boxcar.artery, 'o', MarkerSize=12); hold on;
            plot_over_figure(sig, '--', LineWidth=2);
            plot_over_figure(idl, '-', LineWidth=1.5); hold off;
            ylabel("activity (Bq/mL)")
            xlabel("time (s)")
            saveFigures(pwd, prefix=stackstr(use_dashes=true), closeFigure=true)

            popd(pwd0);

            % deconvBayes
            % boxcar = mlsiemens.BoxcarModel.create( ...
            %     artery=idif_ic, ...
            %     tracer=this.tracer, ...
            %     model_kind=this.model_kind); 
            % boxcar.build_solution();
            % idif_ic = boxcar.artery;
            % idif_ic.fileprefix = mlpipeline.Bids.adjust_fileprefix( ...
            %     idif_ic.fileprefix, ...
            %     post_proc="deconv", ...
            %     remove_substring="_finite");
            % idif_ic.save();
            %figure; plot(idif_ic)
        end
        
        function build_pet_objects(this, opts)
            %% Build PET mips for use as guides when calling draw_tof_mips.
            %  Also build transformations for registering centerlines to native PET.
            %  opts.search_180 == true ~ searchr{x,y,z} = [-180,180]; useful for head-first CO from e7

            arguments 
                this mlaif.MipIdif
                opts.search_180 logical = false
            end

            %% Register T1w on PET avgt.

            if opts.search_180
                t1w_on_pet_flirt = mlfsl.Flirt( ...
                    'in', this.t1w, ...
                    'ref', this.pet_avgt, ...
                    'out', this.t1w_on_pet_fqfn, ...
                    'omat', this.mat(this.t1w_on_pet_fqfn), ...
                    'bins', 256, ...
                    'cost', 'mutualinfo', ...
                    'dof', 6, ...
                    'searchrx', [-180,180], ...
                    'searchry', [-180,180], ...
                    'searchrz', [-180,180], ...
                    'interp', 'spline', ...
                    'noclobber', false);
            else
                t1w_on_pet_flirt = mlfsl.Flirt( ...
                    'in', this.t1w, ...
                    'ref', this.pet_avgt, ...
                    'out', this.t1w_on_pet_fqfn, ...
                    'omat', this.mat(this.t1w_on_pet_fqfn), ...
                    'bins', 256, ...
                    'cost', 'mutualinfo', ...
                    'dof', 6, ...
                    'interp', 'spline', ...
                    'noclobber', false);
            end
            if ~isfile(this.t1w_on_pet_fqfn)
                % do expensive coreg.
                t1w_on_pet_flirt.flirt();
            end
            assert(isfile(this.t1w_on_pet_fqfn))

            %% Transform tof to PET avg, if tof available.

            try
                if ~isempty(this.tof) && isfile(this.tof.fqfn)
                    tof_on_pet_avgt_flirt = copy(t1w_on_pet_flirt);
                    tof_on_pet_avgt_flirt.concatXfm(AtoB=this.mat(this.tof_on_t1w));
                    tof_on_pet_avgt_flirt.in = this.tof;
                    tof_on_pet_avgt_flirt.ref = this.pet_avgt;
                    tof_on_pet_avgt_flirt.out = this.tof_on_pet_fqfn;
                    if ~isfile(this.tof_on_pet_fqfn)
                        tof_on_pet_avgt_flirt.applyXfm();
                    end
                    assert(isfile(this.tof_on_pet_fqfn))
                end
            catch ME
                handwarning(ME)
            end

            this.pet_mipt; % saves
        end
        
        function idif_ic = build_post_centerline(this, opts)
            %% Args:
            %  opts.delete_large_files logical = false
            %  opts.steps logical = true(1, 4)
            %
            %       step(1) ~ calls back_project()
            %       step(2) ~ calls align_centerline_on_tracer()
            %       step(3) ~ calls build_aif(this.pet_dyn)
            %       step(4) ~ calls build_deconv()

            arguments
                this mlaif.MipIdif
                opts.steps logical = true(1, 4)
                opts.delete_large_files logical = false
                opts.frac_select double = this.ALPHA
            end

            % cached on filesystem
            if isfile(this.new_fqfp + ".nii.gz")
                idif_ic = mlfourd.ImagingContext2(this.new_fqfp + ".nii.gz");
                return
            end

            idif_ic = [];
            if opts.steps(1)
                this.back_project(); % overwrites centerline_on_pet.nii.gz
            end
            if opts.steps(2)
                if contains(this.tracer, "co", IgnoreCase=true) || ...
                        contains(this.tracer, "oc", IgnoreCase=true)
                    this.align_centerline_on_tracer(frames=1:25);
                else
                    this.align_centerline_on_tracer();
                end
            end
            if opts.steps(3)
                idif_ic = this.build_aif(this.pet_dyn, frac_select=opts.frac_select); % by statistical sampling
            end
            if opts.steps(4) % && ~strcmpi(this.tracer, "co") && ~strcmpi(this.tracer, "oc")
                idif_ic = this.build_deconv(idif_ic);
            end

            if opts.delete_large_files && ~contains(this.pet_dyn.filepath, "sourcedata")
                deleteExisting(this.pet_dyn.fqfp+".*")
            end
        end

        function idif_ic = build_centerline(this, opts)
            %% Args:
            %  opts.delete_large_files logical = false
            %  opts.steps logical = true

            arguments
                this mlaif.MipIdif
                opts.steps logical = true
                opts.delete_large_files logical = false
                opts.frac_select double = this.ALPHA
            end

            % cached on filesystem
            if isfile(this.new_fqfp + ".nii.gz")
                idif_ic = mlfourd.ImagingContext2(this.new_fqfp + ".nii.gz");
                return
            end

            idif_ic = [];
            if any(opts.steps)
                if isfile(this.tof_on_pet.fqfn)
                    this.draw_tof_centerline();
                else
                    this.draw_t1w_centerline();
                end
            end
        end

        function idif_ic = build_pre_centerline(this, opts)
            %% Args:
            %  opts.delete_large_files logical = false
            %  opts.steps logical = true

            arguments
                this mlaif.MipIdif
                opts.steps logical = true
                opts.delete_large_files logical = false
                opts.frac_select double = this.ALPHA
            end

            % cached on filesystem
            if isfile(this.new_fqfp + ".nii.gz")
                idif_ic = mlfourd.ImagingContext2(this.new_fqfp + ".nii.gz");
                return
            end

            idif_ic = [];
            if any(opts.steps)
                this.build_pet_objects(); % lots of fsl ops
            end
        end
        
        function draw_t1w_mips(this)
            
            pet_miptx = max(this.pet_mipt, [], 1);
            pet_mipty = max(this.pet_mipt, [], 2);
            pet_miptz = max(this.pet_mipt, [], 3);

            pet_miptx.save();
            pet_mipty.save();
            pet_miptz.save();

            %% draw in fsleyes 2-3 voxel widths for barycentric continuity

            pwd0 = pushd(this.pet_avgt.filepath);

            fprintf("Please draw arterial centerlines as overlaid images on the T1w\n" + ...
                "Save files: centerline_zl and centerline_zr.\n\n");
            t1w_mipz = max(this.t1w_on_pet, [], 3); 
            t1w_mipz.save()
            %t1w_mipz.view(pet_miptz) % pet_miptz.view(t1w_mipz)
            fprintf("Please draw arterial centerlines as overlaid images on the T1w\n" + ...
                "Save files: centerline_yl and centerline_yr.\n\n");
            t1w_mipy = max(this.t1w_on_pet, [], 2); 
            t1w_mipy.save()
            %t1w_mipy.view(pet_mipty) % pet_mipty.view(t1w_mipz)
            fprintf("Please draw arterial centerlines as overlaid images on the T1w\n" + ...
                "Save files: centerline_xl and centerline_xr.\n\n");
            t1w_mipx = max(this.t1w_on_pet, [], 1); 
            t1w_mipx.save()
            %t1w_mipx.view(pet_miptx) % pet_miptx.view(t1w_mipz)

            popd(pwd0);
        end
        
        function draw_tof_centerline(this)
            fprintf("Please draw arterial centerlines as overlaid images\n" + ...
                "Save file: centerline_on_pet_nifti2.\n\n");
            this.tof_on_pet.view(this.t1w_on_pet, this.pet_mipt)
        end
        
        function draw_t1w_centerline(this)
            fprintf("Please draw arterial centerlines as overlaid images\n" + ...
                "Save file: centerline_on_pet_nifti2.\n\n");
            this.t1w_on_pet.view(this.pet_mipt)
        end
        
        function draw_tof_mips(this)
            
            pet_miptx = max(this.pet_mipt, [], 1);
            pet_mipty = max(this.pet_mipt, [], 2);
            pet_miptz = max(this.pet_mipt, [], 3);

            pet_miptx.save();
            pet_mipty.save();
            pet_miptz.save();

            %% draw in fsleyes 2-3 voxel widths for barycentric continuity

            pwd0 = pushd(this.pet_avgt.filepath);

            fprintf("Please draw arterial centerlines as overlaid images on the tof\n" + ...
                "Save files: centerline_zl and centerline_zr.\n\n");
            tof_mipz = max(this.tof_on_pet, [], 3); 
            tof_mipz.save()
            %tof_mipz.view(pet_miptz) % pet_miptz.view(tof_mipz)
            fprintf("Please draw arterial centerlines as overlaid images on the tof\n" + ...
                "Save files: centerline_yl and centerline_yr.\n\n");
            tof_mipy = max(this.tof_on_pet, [], 2); 
            tof_mipy.save()
            %tof_mipy.view(pet_mipty) % pet_mipty.view(tof_mipy)
            fprintf("Please draw arterial centerlines as overlaid images on the tof\n" + ...
                "Save files: centerline_xl and centerline_xr.\n\n");
            tof_mipx = max(this.tof_on_pet, [], 1); 
            tof_mipx.save()
            %tof_mipx.view(pet_miptx) % pet_miptx.view(tof_mipx)

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
            source = this.pet_avgt;
            source_on_target = mlfourd.ImagingContext2( ...
                source.fqfp+"_on_"+target_trc+".nii.gz");
            cl_flirt = mlfsl.Flirt( ...
                'in', source, ...
                'ref', target, ...
                'out', source_on_target, ...
                'omat', this.mat(source_on_target), ...
                'bins', 256, ...
                'cost', this.FLIRT_COST, ...
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
        
        function ic = make_volume(this, img, opts)
            arguments
                this mlaif.MipIdif %#ok<INUSA>
                img {mustBeNumeric}
                opts.ic_template mlfourd.ImagingContext2 = this.pet_avgt
                opts.fileprefix {mustBeTextScalar} = stackstr(3)
            end
            ifc = opts.ic_template.nifti;
            ifc.img = img;
            ifc.fileprefix = opts.fileprefix;
            ic = mlfourd.ImagingContext2(ifc);
        end  
        
        function g = new_filepath(this)
            g = myfileparts(this.new_fqfp);
        end
        
        function fqfp = new_fqfp(this, opts)
            %% seeks out pertinent fqfp in derivatives, not sourcedata

            arguments
                this mlaif.MipIdif
                opts.remove_substring {mustBeTextScalar} = "_timeAppend-4"
                opts.tags {mustBeTextScalar} = ""
            end

            pth = fullfile(this.mediator_.derivSesPath, "pet");
            fp = mlpipeline.Bids.adjust_fileprefix(this.pet_avgt.fileprefix, ...                
                new_proc="MipIdif"+opts.tags, new_mode="idif", remove_substring=opts.remove_substring);
            fqfp = fullfile(pth, fp);
            if contains(fqfp, "*")
                fqfp = mglob(fqfp);
                assert(~isempty(fqfp))
                if numel(fqfp) > 1
                    warning("mlaif:RuntimeWarning", stackstr()+" returned string array of length "+numel(fqfp))
                end
            end
        end
        
        function g = new_fqfileprefix(this, varargin)
            g = this.new_fqfp(varargin{:});
        end
        
        function pet_dyn_on_anatomy(~)
            error("mlaif:NotImplementedError", "MipIdif");
        end
        
        function idif_ic = save_imaging_context(this, idif, opts)
            %% adjusts fileprefix and json_metadata prior to saving

            arguments
                this mlaif.MipIdif
                idif mlfourd.ImagingContext2
                opts.filepath = this.pet_avgt.filepath
                opts.fileprefix_template = this.pet_avgt.fileprefix
                opts.json_metadata_template = this.pet_dyn.json_metadata
                opts.remove_substring {mustBeTextScalar} = "_timeAppend-4"
                opts.tags {mustBeTextScalar} = ""
            end
            if opts.tags ~= "" && ~strcmp(opts.tags{1}(1), "-")
                opts.tags = "-"+opts.tags;
            end

            idif_ic = mlfourd.ImagingContext2(idif);
            idif_ic.filepath = opts.filepath;

            % adjust fileprefix from template
            fp = opts.fileprefix_template;
            idif_ic.fileprefix = mlpipeline.Bids.adjust_fileprefix(fp, ...                
                new_proc="MipIdif"+opts.tags, new_mode="idif", remove_substring=opts.remove_substring);

            % adjust json_metadata from template
            j = opts.json_metadata_template;
            j = mlpipeline.ImagingMediator.ensureNumericTimingData(j);
            idif_ic.json_metadata = j;
            idif_ic.selectNiftiTool();
            idif_ic.save();
            this.product_ = idif_ic;
        end
        
        function ic = select_frames(this, dyn_ic, opts)
            arguments
                this mlaif.MipIdif
                dyn_ic mlfourd.ImagingContext2
                opts.frames double {mustBePositive}
                opts.stat function_handle = @max
            end
            if isscalar(opts.frames)
                opts.frames = 1:opts.frames;
            end
            ifc = dyn_ic.imagingFormat;
            ifc.img = ifc.img(:,:,:,opts.frames);
            try
                ifc.img = opts.stat(ifc.img, [], 4);
            catch ME
                fprintf("%s:  %s\n", stackstr(), ME.message);
                ifc.img = opts.stat(ifc.img, 4);
            end
            ifc.fileprefix = sprintf("%s_frames%i-%i_%s", ...
                ifc.fileprefix, opts.frames(1), opts.frames(end), func2str(opts.stat));
            ic = mlfourd.ImagingContext2(ifc);
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
            aife.fileprefix = sprintf("aif_trc-%s_proc-%s", trc, stackstr(use_dashes=true));
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
            aifl.fileprefix = sprintf("aif_trc-%s_proc-%se", trc, stackstr(use_dashes=true));
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
        function ic = apply_mask(ic_in, ic_msk, blur, opts)
            arguments
                ic_in mlfourd.ImagingContext2
                ic_msk mlfourd.ImagingContext2
                blur {mustBeScalarOrEmpty} = 5  % 5 mm fwhh blur; flirt fails with 40 mm
                opts.binarize logical = false
            end

            fileprefix0 = ic_in.fileprefix;
            ifc = ic_in.imagingFormat;
            if ~isempty(blur)
                ic_msk1 = ic_msk.blurred(blur);
            else
                ic_msk1 = ic_msk;
            end
            if opts.binarize
                ic_msk1 = ic_msk1.thresh(0.1);
                ic_msk1 = ic_msk1.binarized();
            end
            ifc.img = ifc.img .* ic_msk1.imagingFormat.img;
            
            % package product
            ic = mlfourd.ImagingContext2(ifc);
            ic.fileprefix = fileprefix0 + "_apply-mask";
            ic.save();
        end
        
        function this = create(opts)
            %% bids_kit has file in sourceSesPet

             arguments
                opts.bids_kit mlkinetics.BidsKit
                opts.tracer_kit mlkinetics.TracerKit
                opts.scanner_kit mlkinetics.ScannerKit
                opts.pet_avgt = []
                opts.pet_mipt = []
                opts.model_kind = ""
                opts.reference_tracer {mustBeTextScalar} = "ho" % consider "ho" or "fdg"
                opts.minz_for_mip = 5 % used to remove scattering from inferior FOV
                opts.dilate_m double = 0
                opts.filename_pattern_dynamic {mustBeTextScalar} = "sub-*_ses-*delay*BrainMoCo2*createNiftiMovingAvgFrames.nii.gz"
                opts.filename_pattern_static {mustBeTextScalar} = "sub-*_ses-*delay*BrainMoCo2*createNiftiStatic.nii.gz"
                opts.filename_pattern_dynamic_fdg {mustBeTextScalar} = "sub-*_ses-*_trc-fdg_proc-consoleDynamic.nii.gz"
                opts.filename_pattern_static_fdg {mustBeTextScalar} = "sub-*_ses-*_trc-fdg_proc-consoleStatic.nii.gz"
            end

            this = mlaif.MipIdif();
            this.bids_kit_ = opts.bids_kit;
            this.tracer_kit_ = opts.tracer_kit;
            this.scanner_kit_ = opts.scanner_kit;
            
            this.mediator_ = this.bids_kit_.make_bids_med();
            this.t1w_ = this.mediator_.t1w_ic;
            try
                this.tof_ = this.mediator_.tof_ic;
            catch ME
                fprintf("%s: tof not found\n", stackstr());
            end

            % pet_dyn:
            % defer to filesystem this expensive object;
            % avoid storing in memory because large size will force
            % ImagingContext2 to return handles rather than safe copies
            assert(isfile(this.mediator_.fqfn))
            assert(contains(this.mediator_.filepath, this.mediator_.sourceSesPath))

            if ~isempty(opts.pet_avgt)
                this.pet_avgt_ = mlfourd.ImagingContext2(opts.pet_avgt);
            end
            if ~isempty(opts.pet_mipt)
                this.pet_mipt_ = mlfourd.ImagingContext2(opts.pet_mipt);
            end
            this.model_kind_ = opts.model_kind;
            if isemptytext(opts.model_kind)
                if strcmpi(this.tracer, "OO")
                    this.model_kind_ = "4bolus";
                else
                    this.model_kind_ = "3bolus";
                end
            end
            this.reference_tracer = opts.reference_tracer;
            this.minz_for_mip = opts.minz_for_mip;
            this.dilate_m_ = opts.dilate_m;
            this.filename_pattern_dynamic = opts.filename_pattern_dynamic;
            this.filename_pattern_static = opts.filename_pattern_static;
            this.filename_pattern_dynamic_fdg = opts.filename_pattern_dynamic_fdg;
            this.filename_pattern_static_fdg = opts.filename_pattern_static_fdg;
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
        
        function ic = timeAppend(ic)
            toglob = strrep(ic, "-delay0-", "-delay*-");
            mg = mglob(toglob);
            if 1 == length(mg)
                return
            else
                fprintf("%s: (%s).timeAppend(%s)\n", stackstr(3), ic.fileprefix, mg(2))
                ic = ic.timeAppend(mg(2));
            end
        end
        
        function viewAnat()
        end
    end

    %% PROTECTED

    properties (Access = protected)
        centerline_on_reftrc_
        centerline_on_pet_
        dilate_m_
        fdg_avgt_
        fdg_dyn_
        json_metadata_template_
        mediator_
        model_kind_
        pet_avgt_
        pet_dyn_
        pet_mipt_
        product_
        reftrc_dyn_
        reftrc_avgt_
        t1w_
        t1w_on_pet_
        t1w_on_pet_fqfn_
        tof_
        tof_on_pet_
        tof_on_pet_fqfn_
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

            if ~isempty(this.fdg_avgt_)
                that.fdg_avgt_ = copy(this.fdg_avgt_); end
            if ~isempty(this.fdg_dyn_)
                that.fdg_dyn_ = copy(this.fdg_dyn_); end

            if ~isempty(this.pet_avgt_)
                that.pet_avgt_ = copy(this.pet_avgt_); end
            if ~isempty(this.pet_dyn_)
                that.pet_dyn_ = copy(this.pet_dyn_); end
            if ~isempty(this.pet_mipt_)
                that.pet_mipt_ = copy(this.pet_mipt_); end

            if ~isempty(this.reftrc_dyn_)
                that.reftrc_dyn_ = copy(this.reftrc_dyn_); end
            if ~isempty(this.reftrc_avgt_)
                that.reftrc_avgt_ = copy(this.reftrc_avgt_); end

            if ~isempty(this.t1w_)
                that.t1w_ = copy(this.t1w_); end
            if ~isempty(this.t1w_on_pet_)
                that.t1w_on_pet_ = copy(this.t1w_on_pet_); end

            if ~isempty(this.tof_)
                that.tof_ = copy(this.tof_); end
            if ~isempty(this.tof_on_pet_)
                that.tof_on_pet_ = copy(this.tof_on_pet_); end

            if ~isempty(this.bids_kit_)
                that.bids_ = copy(this.bids_kit_); end            
            if ~isempty(this.tracer_kit_)
                that.bids_ = copy(this.tracer_kit_); end
            if ~isempty(this.scanner_kit_)
                that.bids_ = copy(this.scanner_kit_); end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
