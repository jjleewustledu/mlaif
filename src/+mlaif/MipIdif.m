classdef MipIdif < handle & mlsystem.IHandle
    %% is a builder.  Start with build_all(this).  
    %  
    %  Created 12-Jun-2023 22:59:36 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlaif/src/+mlaif.
    %  Developed on Matlab 9.14.0.2254940 (R2023a) Update 2 for MACI64.  Copyright 2023 John J. Lee.
    
    properties (Constant)
        ALPHA = 0.1
        max_len_mipt = 30
    end

    properties
        minz_for_mip
        reference_tracer
    end
    
    properties (Dependent)
        centerline_on_reftrc
        centerline_on_pet
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
        tof
        tof_on_pet
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
                        fullfile(pth, "centerline_on_pet.nii.gz"));
                catch ME
                    handexcept(ME)
                end
            end
            g = copy(this.centerline_on_reftrc_);
        end
        function g = get.centerline_on_pet(this)
            if isempty(this.centerline_on_pet_)
                try
                    this.centerline_on_pet_ = mlfourd.ImagingContext2( ...
                        fullfile(this.pet_avgt.filepath, "centerline_on_pet.nii.gz"));
                catch ME
                    handexcept(ME)
                end
            end
            g = copy(this.centerline_on_pet_);
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

            % prefer static
            pth = fullfile(this.mediator_.sourceSesPath, "pet");
            mg = mglob(fullfile(pth, "sub-*_ses-*delay*BrainMoCo2*createNiftiStatic.nii.gz"));
            if ~isemptytext(mg) && isfile(mg(1))
                % pick delay0
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
            this.pet_avgt_ = this.pet_dyn.timeAveraged(weights=this.pet_dyn.taus);
            this.pet_avgt_.relocateToDerivativesFolder();
            if ~isfile(this.pet_avgt_)
                save(this.pet_avgt_);
            end
            g = this.pet_avgt_;
        end
        function g = get.pet_dyn(this)
            %% pet_dyn will commonly be too large for in-memory copying; 
            %  therefor, always instantiate as FilesystemTools in ImagingContext2.  

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
            L = min(size(ifc, 4), this.max_len_mipt);
            ifc.img = ifc.img(:,:,this.minz_for_mip:end,1:L);
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
            mg = mglob(fullfile(this.mediator_.sourceSubPath, "ses-*", "pet", "sub-*_ses-*_trc-"+this.reference_tracer+"_proc-delay*BrainMoCo2-createNiftiMovingAvgFrames.nii.gz"));
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
            mg = mglob(fullfile(this.mediator_.sourceSubPath, "ses-*", "pet", "sub-*_ses-*_trc-"+this.reference_tracer+"_proc-delay*BrainMoCo2-createNiftiStatic.nii.gz"));
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
        function g = get.tof(this)
            g = this.mediator_.tof_ic;
        end
        function g = get.tof_on_pet(this)
            if ~isempty(this.tof_on_pet_)
                g = this.tof_on_pet_;
                return
            end

            if isfile(this.tof_on_pet_fqfn_)
                this.tof_on_pet_ = mlfourd.ImagingContext2(this.tof_on_pet_fqfn_);
                g = this.tof_on_pet_;
                return
            end

            this.build_pet_objects();
            g = mlfourd.ImagingContext2(this.tof_on_pet_fqfn_);
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
            %% aligns FDG centerline on other tracers

            arguments
                this mlaif.MipIdif
                opts.frames double = []
            end

            if contains(this.pet_dyn_fqfn_, "_trc-"+this.reference_tracer, IgnoreCase=true) % ho is in "/home"
                cl = this.centerline_on_pet;
                assert(isfile(cl));
                return
            end

            % use early ref frames to align vasculature

            if ~isempty(opts.frames)
                reftrac_avgt__ = this.select_frames(this.reftrc_dyn, frames=opts.frames);
            else
                reftrac_avgt__ = this.reftrc_avgt;
            end

            % use dlicv masks

            inweight = this.build_dlicv_mask(this.pet_avgt);
            refweight = this.build_dlicv_mask(reftrac_avgt__);

            pet_on_ref_fqfn = this.pet_avgt.fqfp+"_on_"+this.reference_tracer+".nii.gz"; 
            pet_on_ref = mlfourd.ImagingContext2(pet_on_ref_fqfn);
            pet_on_ref_flirt = mlfsl.Flirt( ...
                'in', this.pet_avgt, ...
                'inweight', inweight, ...
                'ref', reftrac_avgt__, ...
                'refweight', refweight, ...
                'out', pet_on_ref, ...
                'omat', this.mat(pet_on_ref), ...
                'bins', 256, ...
                'cost', 'mutualinfo', ...
                'dof', 6, ...
                'interp', 'spline', ...
                'noclobber', false);
            if ~isfile(pet_on_ref.fqfn)
                % do expensive coreg.
                pet_on_ref_flirt.flirt();
            end
            assert(isfile(pet_on_ref.fqfn))

            pet_on_ref_flirt.invertXfm();

            pet_on_ref_flirt.in = this.centerline_on_reftrc;
            pet_on_ref_flirt.out = this.centerline_on_pet;
            pet_on_ref_flirt.interp = "trilinear";
            pet_on_ref_flirt.applyXfm();
            copyfile(this.centerline_on_pet.fqfn, this.centerline_on_pet.fqfp+"_trilinear.nii.gz")
            this.centerline_on_pet_ = this.centerline_on_pet.thresh(this.ALPHA);
            this.centerline_on_pet_ = this.centerline_on_pet.binarized();
            this.centerline_on_pet_.fileprefix = "centerline_on_pet";
            this.centerline_on_pet_.save();
            assert(isfile(this.centerline_on_pet))
        end
        function ic = build_dlicv_mask(this, pet_static)
            arguments
                this mlaif.MipIdif
                pet_static {mustBeNonempty}
            end
            pet_static = mlfourd.ImagingContext2(pet_static);

            t1w_on_pet_fqfn = fullfile(pet_static.filepath, "T1w_on_"+pet_static.filename); 
            t1w_on_pet = mlfourd.ImagingContext2(t1w_on_pet_fqfn);
            t1w_on_pet_flirt = mlfsl.Flirt( ...
                'in', this.t1w, ...
                'ref', pet_static, ...
                'out', t1w_on_pet, ...
                'omat', this.mat(t1w_on_pet), ...
                'bins', 256, ...
                'cost', 'mutualinfo', ...
                'dof', 6, ...
                'interp', 'spline', ...
                'noclobber', false);
            if ~isfile(t1w_on_pet.fqfn) || ~isfile(this.mat(t1w_on_pet))
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
            cl.fileprefix = 'centerline_on_pet';
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
            
            % adjust opts.mipt_thr for CO
            tr = mlfourd.ImagingContext2(petobj);
            re = mlpipeline.Bids.regexp_fileprefix(tr);
            assert(~isempty(re))
            if strcmp(re.trc, "co") || strcmp(re.trc, "oc")
                opts.mipt_thr = min(50000, opts.mipt_thr);
            end

            %% select 0.05*numel of brightest and earliest-arriving voxels from centerline        

            centerline = logical(opts.cl);
            cache.centerline = centerline;

            % tr_mipt, tr_indices ~ MIP, time of max
            tr_single = single(tr);
            sz = size(tr_single);
            [tr_mipt,tr_indices] = max(tr, [], 4);
            tr_mipt.save();
            tr_indices.save();
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
            T1 = T(ascol(1:round(opts.frac_select*len)), :);
            cache.T1 = T1;
            tr_single___ = tr_single__(T1.indices,:);
            cache.tr_single___ = tr_single___;

            % save cache
            cache_fqfn = sprintf("%s_%s_cache.mat", this.new_fqfp, stackstr());
            save(cache_fqfn, "cache");
            
            % aif from volume-average of brightest 0.05 from centerline
            idif_ic = this.save_imaging_context(mlfourd.ImagingContext2(mean(tr_single___, 1)));

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
            arguments
                this mlaif.MipIdif
                opts.steps logical = true(1, 6)
                opts.delete_large_files logical = false;
            end

            % cached on filesystem
            if isfile(this.new_fqfp + ".nii.gz")
                idif_ic = mlfourd.ImagingContext2(this.new_fqfp + ".nii.gz");
                return
            end

            idif_ic = [];
            if ~isfile(this.centerline_on_pet.fqfn)
                if opts.steps(1)
                    this.build_pet_objects(); % lots of fsl ops
                end
                if opts.steps(2)
                    this.build_tof_mips(); % write intermediates needed for drawing on MIPs
                end
                if opts.steps(3)
                    this.back_project(); % overwrites centerline_on_pet.nii.gz
                end
            end
            if opts.steps(4)
                if contains(this.tracer, "co", IgnoreCase=true) || ...
                        contains(this.tracer, "oc", IgnoreCase=true)
                    this.align_centerline_on_tracer(frames=1:25);
                else
                    this.align_centerline_on_tracer();
                end
            end
            if opts.steps(5)
                this.pet_dyn_fqfn_ = mlkinetics.BidsKit.timeAppendForFqfn( ...
                    this.pet_dyn_fqfn_, ...
                    do_save=true);
                idif_ic = this.build_aif(this.pet_dyn); % by statistical sampling
            end
            if opts.steps(6) % && ~strcmpi(this.tracer, "co") && ~strcmpi(this.tracer, "oc")
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
        function build_pet_objects(this)
            %% Build PET mips for use as guides when calling build_tof_mips.
            %  Also build transformations for registering centerlines to native PET.

            %% Register T1w on PET avgt.

            t1w_on_pet_fqfn = fullfile(this.pet_avgt.filepath, "T1w_on_"+this.pet_avgt.fileprefix+".nii.gz"); 
            t1w_on_pet = mlfourd.ImagingContext2(t1w_on_pet_fqfn);
            t1w_on_pet_flirt = mlfsl.Flirt( ...
                'in', this.t1w, ...
                'ref', this.pet_avgt, ...
                'out', t1w_on_pet, ...
                'omat', this.mat(t1w_on_pet), ...
                'bins', 256, ...
                'cost', 'mutualinfo', ...
                'dof', 6, ...
                'interp', 'spline', ...
                'noclobber', false);
            if ~isfile(t1w_on_pet.fqfn)
                % do expensive coreg.
                t1w_on_pet_flirt.flirt();
            end
            assert(isfile(t1w_on_pet.fqfn))

            %% Transform tof to PET avgt.

            this.tof_on_pet_ = mlfourd.ImagingContext2(this.tof_on_pet_fqfn_);
            tof_on_pet_avgt_flirt = copy(t1w_on_pet_flirt);
            tof_on_pet_avgt_flirt.concatXfm(AtoB=this.mat(this.tof_on_t1w));
            tof_on_pet_avgt_flirt.in = this.tof;
            tof_on_pet_avgt_flirt.ref = this.pet_avgt;
            tof_on_pet_avgt_flirt.out = this.tof_on_pet_;
            if ~isfile(this.tof_on_pet_.fqfn)
                tof_on_pet_avgt_flirt.applyXfm();
            end
            assert(isfile(this.tof_on_pet_.fqfn))

            this.pet_mipt; % saves
        end
        function build_tof_mips(this)
            
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
            end

            pth = fullfile(this.mediator_.derivSesPath, "pet");
            fp = mlpipeline.Bids.adjust_fileprefix(this.pet_avgt.fileprefix, ...                
                new_proc="MipIdif-finite-deconv", new_mode="idif", remove_substring=opts.remove_substring);
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
            end

            idif_ic = mlfourd.ImagingContext2(idif);
            idif_ic.filepath = opts.filepath;

            % adjust fileprefix from template
            fp = opts.fileprefix_template;
            idif_ic.fileprefix = mlpipeline.Bids.adjust_fileprefix(fp, ...                
                new_proc="MipIdif", new_mode="idif", remove_substring="_timeAppend-4");

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
            end
            if isscalar(opts.frames)
                opts.frames = 1:opts.frames;
            end
            ifc = dyn_ic.imagingFormat;
            ifc.img = ifc.img(:,:,:,opts.frames);
            ifc.img = mean(ifc.img, 4);
            ifc.fileprefix = sprintf("%s_frames%i-%i", ifc.fileprefix, opts.frames(1), opts.frames(end));
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
        function this = create(opts)
            %% bids_kit has file in sourceSesPet

             arguments
                opts.bids_kit mlkinetics.BidsKit
                opts.tracer_kit mlkinetics.TracerKit
                opts.scanner_kit mlkinetics.ScannerKit
                opts.pet_avgt = []
                opts.pet_mipt = []
                opts.model_kind = ""
                opts.reference_tracer {mustBeTextScalar} = "fdg" % consider "ho"
                opts.minz_for_mip = 10 % used to remove scattering from inferior FOV
            end

            this = mlaif.MipIdif();
            this.bids_kit_ = opts.bids_kit;
            this.tracer_kit_ = opts.tracer_kit;
            this.scanner_kit_ = opts.scanner_kit;
            
            this.mediator_ = this.bids_kit_.make_bids_med();
            this.t1w_ = this.mediator_.t1w_ic;
            this.tof_ = this.mediator_.tof_ic;

            % pet_dyn:
            % defer to filesystem this expensive object;
            % avoid storing in memory because large size will force
            % ImagingContext2 to return handles rather than safe copies
            assert(isfile(this.mediator_.fqfn))
            assert(contains(this.mediator_.filepath, this.mediator_.sourceSesPath))
            this.pet_dyn_fqfn_ = this.mediator_.fqfn;

            % tof_on_pet:  
            % defer to filesystem this expensive object
            this.tof_on_pet_fqfn_ = fullfile(this.pet_avgt.filepath, "tof_on_"+this.pet_avgt.fileprefix+".nii.gz"); 

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
        json_metadata_template_
        mediator_
        model_kind_
        pet_avgt_
        pet_dyn_fqfn_
        pet_mipt_
        product_
        reftrc_dyn_
        reftrc_avgt_
        t1w_
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
            if ~isempty(this.pet_avgt_)
                that.pet_avgt_ = copy(this.pet_avgt_); end
            if ~isempty(this.pet_mipt_)
                that.pet_mipt_ = copy(this.pet_mipt_); end
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
