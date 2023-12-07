classdef Idif < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 26-Sep-2023 15:56:00 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlaif/src/+mlaif.
    %  Developed on Matlab 9.14.0.2337262 (R2023a) Update 5 for MACI64.  Copyright 2023 John J. Lee.
    
    properties   
    end

    properties (Dependent)
        Measurement
 		timeInterpolants    
        tracer
    end

    methods %% GET
        function g = get.Measurement(this)
            g = this.Measurement_;
        end
        function g = get.timeInterpolants(this)
            if ~isempty(this.timeInterpolants_)
                g = this.timeInterpolants_;
                return
            end

            try
                dt = this.json_.starts(2) - this.json_.starts(1);
                timesF = this.json_.timesMid(end);
                this.timeInterpolants_ = 0:dt:timesF(end);
                g = this.timeInterpolants_;
            catch ME
                handexcept(ME)
            end
        end
        function g = get.tracer(this)
            g = this.tracer_;
        end
    end

    methods
        function q = deconv(this, opts)
            arguments
                this mlaif.Idif
                opts.return_ic logical = true
            end

            ifc__ = this.Measurement.imagingFormat;
            L = this.moving_aver_oper(tracer=this.tracer);
            q = asrow(lsqnonneg(L, ascol(double(ifc__.img))));
            q = q(1:size(L,1));
            ifc__.img = q;
            ifc__.fileprefix = ifc__.fileprefix+"-deconv";
            ifc__.save();

            if opts.return_ic
                q = mlfourd.ImagingContext2(ifc__);
            end            
        end
        function q = deconvBayes(this, opts)
            arguments
                this mlaif.Idif
                opts.return_ic logical = true
            end

            ifc__ = this.Measurement.imagingFormat;
            ral = mlswisstrace.RadialArteryLee2021( ...
                'tracer', convertStringsToChars(this.tracer), ...
                'kernel', this.kernel, ...
                'model_kind', '3bolus', ...
                'Measurement', asrow(double(this.Measurement)));
            ral = ral.solve();
            ral.plot()
            q = ral.deconvolved();
            ifc__.img = q;
            ifc__.fileprefix = ifc__.fileprefix+"-deconvBayes";
            ifc__.save();

            if opts.return_ic
                q = mlfourd.ImagingContext2(ifc__);
            end     
        end
        function q = deconvBayes2(this, opts)
            arguments
                this mlaif.Idif
                opts.return_ic logical = true
            end

            ifc__ = this.Measurement.imagingFormat;
            j = ifc__.json_metadata;
            timesMid = ascol(j.timesMid);
            activityDensity = ascol(ifc__.img);
            T = table(timesMid, activityDensity);
            art = mlaif.Artery( ...
                T, ...
                'tracer', convertStringsToChars(this.tracer), ...
                'kernel', this.kernel, ...
                'model_kind', '3bolus');
            art = art.solve();
            art.plot()
            q = art.deconvolved();
            ifc__.img = q;
            ifc__.fileprefix = ifc__.fileprefix+"-deconvBayes2";
            ifc__.save();

            if opts.return_ic
                q = mlfourd.ImagingContext2(ifc__);
            end     
        end
        function q = deconvBayes3(this, opts)
            arguments
                this mlaif.Idif
                opts.return_ic logical = true
            end

            ifc__ = this.Measurement.imagingFormat;
            j = ifc__.json_metadata;
            timesMid = ascol(j.timesMid);
            activityDensity = ascol(ifc__.img);
            T = table(timesMid, activityDensity);
            art = mlaif.Artery( ...
                T, ...
                'tracer', convertStringsToChars(this.tracer), ...
                'kernel', 1, ...
                'model_kind', '3bolus-window10');
            art = art.solve();
            art.plot()
            q = art.deconvolved();
            ifc__.img = q;
            ifc__.fileprefix = ifc__.fileprefix+"-deconvBayes2";
            ifc__.save();

            if opts.return_ic
                q = mlfourd.ImagingContext2(ifc__);
            end     
        end
        function k = kernel(this)
            k = zeros(size(this.timeInterpolants));
            Nk = this.json_.taus(1);
            t0 = 0; % round(this.json_.taus(1)/2);
            k(t0+1:t0+Nk) = 1/Nk;
        end
        function plot_deconv(this)
        end
    end

    methods (Static)
        function idif_ic = centerline2idif(cl, petdyn, opts)
            arguments
                cl {mustBeNonempty}
                petdyn {mustBeNonempty}
                opts.thresh double = 180000
            end
            cl = mlfourd.ImagingContext2(cl);
            petdyn = mlfourd.ImagingContext2(petdyn);
            dest_path = strrep(petdyn.filepath, "souredata", "derivatives");
            ensuredir(dest_path);

            % parse petdyn for tags
            re = mlaif.Idif.regexp_for_tags();
            if strcmp(re.trc, "co") || strcmp(re.trc, "oc")
                opts.thresh = min(50000, opts.thresh);
            end

            % idif from volume-averaged centerline
            idif_ic = petdyn.volumeAveraged(cl);
            idif_ic.fileprefix = sprintf("%s_%s_trc-%s_proc-%s", re.sub, re.ses, re.trc, stackstr(use_dashes=true)) ;
            idif_ic.filepath = dest_path;
            idif_ic.save();

            % idif deconvolved from moving averages
            L = mlaif.Idif.moving_aver_oper(tracer=re.trc);
            idif_ifc = idif_ic.imagingFormat;
            img__ = asrow(lsqnonneg(L, ascol(double(idif_ifc.img))));
            idif_ifc.img = img__(1:size(L,1));
            idif_ifc.fileprefix = idif_ifc.fileprefix+"-lsqnoneg-Lma";
            idif_ifc.save();

            % plot
            h = figure;
            timesMid = asrow(idif_ic.json_metadata.timesMid);
            plot( ...
                timesMid, asrow(idif_ic.imagingFormat.img), ...
                timesMid, asrow(idif_ifc.img), "o", MarkerSize=8)
            fontsize(scale=1.3)
            legend(["moving-average", "lsqnoneg L_{ma}"])
            xlabel("times (s)")
            ylabel("IDIF activity density (Bq/mL)")
            title(sprintf("IDIF by centerline, %s", idif_ic.fileprefix), Interpreter="none")
            saveFigure2(h, idif_ic.fqfp)
        end
        function this = create(Measurement, opts)
            arguments
                Measurement {mustBeNonempty}
                opts.centerline = []
                opts.json struct = struct([])
                opts.tracer {mustBeText} = ""
            end
            Measurement = mlfourd.ImagingContext2(Measurement);
            re = mlaif.Idif.regexp_for_tags(Measurement);
            if ~isempty(opts.centerline) && ndims(Measurement) == 4
                fprintf("%s: volume-averaging %s with centerline... \n", ...
                    stackstr(), Measurement.fileprefix);
                ic__ = Measurement.volumeAveraged(mlfourd.ImagingContext2(opts.centerline));
                ic__.fileprefix = sprintf("%s_%s_trc-%s_proc-%s", re.sub, re.ses, re.trc, stackstr(use_dashes=true));
                ic__.save();
                Measurement = ic__;
            end

            this = mlaif.Idif();
            this.Measurement_ = Measurement;
            if ~isempty(this.Measurement_.json_metadata)
                this.json_ = this.Measurement_.json_metadata;
            end
            if ~isempty(opts.json)
                this.json_ = opts.json;
            end
            if isfield(re, "trc")
                this.tracer_ = re.trc;
            end
            if ~isemptytext(opts.tracer)
                this.tracer_ = opts.tracer;
            end
        end

        function f1 = moving_average(f, dT)
            arguments
                f double
                dT double = 10
            end
            f = asrow(f);
            Nf = length(f);
            Nf1 = Nf - dT + 1;
            f1 = zeros([1, Nf1]);
            for ti = 1:Nf1
                f1(ti) = mean(f(ti:ti+dT-1));
            end
        end
        function L = moving_aver_oper(opts)
            arguments
                opts.starts double = 0:9
                opts.taus double = []
                opts.tracer = ""
            end
            if isempty(opts.taus)
                switch convertStringsToChars(opts.tracer)
                    case {'co', 'oc'}
                        opts.taus=10*ones(1,29);
                    case {'ho', 'oo'}
                        opts.taus=10*ones(1,11);
                    case 'fdg'
                        opts.taus=10*ones(1,359);
                    otherwise
                        error("mlkinetics:ValueError", "%s: tracer->%s", stackstr(), this.tracer)
                end
            end
            M = length(opts.starts);
            N = length(opts.taus);
            P = M*N;
            P1 = M*N + M - 1;

            tau = opts.taus(1);
            assert(tau == M, stackstr());

            L = tril(ones(P,P1), M-1) - tril(ones(P,P1), -1);
            L = L/M;
        end
        function re = regexp_for_tags(obj)
            %% re.{sub,ses,trc,proc}

            if isa(obj, "mlfourd.ImagingContext2")
                obj = obj.fileprefix;
            end
            [~,fp] = myfileparts(obj);
            re = regexp(fp, "(?<sub>sub-\S+)_(?<ses>ses-\d+)_trc-(?<trc>\S+)_proc-(?<proc>\S+)", "names");
            assert(~isempty(re))
        end
    end

    %% PRIVATE

    properties (Access = private)
        json_
        Measurement_
        timeInterpolants_
        tracer_
    end

    methods (Access = private)
        function this = Idif()
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
