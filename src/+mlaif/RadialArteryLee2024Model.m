classdef RadialArteryLee2024Model < handle & mlsystem.IHandle
    %% line1
    %  line2
    %  
    %  Created 21-Dec-2023 16:10:30 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlaif/src/+mlaif.
    %  Developed on Matlab 23.2.0.2459199 (R2023b) Update 5 for MACA64.  Copyright 2023 John J. Lee.
    
    properties (Constant) 
        ks_names = {'\alpha' '\beta' 'p' 'dp_F' 't_0' 'steadystate\_fraction' 'bolus_2\_fraction' 'bolus_2\_delay' 'recirc\_fraction' 'recirc\_delay' 'amplitude\_fraction' '\gamma'}
    end
    
	properties 	
        closeFigures
        Data        
        LENK
 		map
        Nensemble
        scalingWorldline
        t0_forced
    end

    properties (Dependent)
        artery
        artery_interpolated % double
        kernel
        measurement 
        model
        product % mlfourd.ImagingContext2
        times_sampled % double
    end

    methods %% GET
        function g = get.artery(this)
            g = copy(this.artery_);
        end
        function g = get.artery_interpolated(this)
            if ~isempty(this.artery_interpolated_)
                g = this.artery_interpolated_;
                return
            end

            this.artery_interpolated_ = asrow(this.artery.imagingFormat.img);
            g = this.artery_interpolated_;
        end
        function g = get.kernel(this)
            g = this.Data.kernel;
        end
        function g = get.measurement(this)
            g = this.measurement_;
        end
        function g = get.model(this)
            g = this;
        end
        function g = get.product(this)
            g = copy(this.product_);
        end
        function g = get.times_sampled(this)
            g = this.times_sampled_;
        end
    end

    methods
        function this = build_model(this, opts)
            arguments
                this mlaif.RadialArteryLee2024Model
                opts.map containers.Map = this.preferredMap()
                opts.measurement {mustBeNumeric} = []
                opts.solver_tags = "simulanneal"
            end

            this.map = opts.map;
            this.adjustMapForTracer();
            this.adjustMapForModelKind();
            if ~isempty(opts.measurement)
                this.measurement_ = opts.measurement;
            end
            if contains(opts.solver_tags, "simulanneal")
                this.solver_ = mlswisstrace.RadialArteryLee2021SimulAnneal(context=this);
            end
            if contains(opts.solver_tags, "multinest")
                this.solver_ = mlnest.MultiNest(context=this);
            end
        end
        function soln = build_solution(this)
            %% MAKE_SOLUTION
            %  @return ks_ in R^1 as mlfourd.ImagingContext2, without saving to filesystems.   
            %  @return this.artery, which is deconvolved.

            % this.Data is assigned in create(), since build_solution() requires no adjustments of this.Data

            % solve artery model and insert solutions into ks
            this.build_model(measurement=this.measurement_); % model assigns this.measurement_

            solved = cell(this.Nensemble, 1);
            losses = NaN(this.Nensemble, 1);
            for idx = 1:this.Nensemble
                tic

                % find lowest loss in the ensemble
                solved{idx} = this.solver_solve();
                losses(idx) = solved{idx}.product.loss;

                toc
            end

            % find best loss in solved
            T = table(solved, losses);
            T = sortrows(T, "losses", "ascend");
            solved_star = T{1, "solved"}{1};

            % update idealized artery
            pr_ = solved_star.product;
            [~,A,ideal] = this.sampled(pr_.ks, this.Data, [], solved_star.TimesSampled);
            img_new_ = solved_star.rescaleModelEstimate(this.scalingWorldline*A*ideal);
            ai_new_ = copy(this.artery.selectImagingTool(img=img_new_));
            ai_new_.fileprefix = strrep(this.artery.fileprefix, "_pet", "_radialarterylee");
            this.artery_ = ai_new_;  
            solved_star.ArteryInterpolated = img_new_; %this.interp1_artery(img=img_new_);      
            this.solver_ = solved_star;

            % product_ := ks
            ks_mat_= [asrow(pr_.ks), pr_.loss];
            ks_mat_ = single(ks_mat_);
            soln = this.artery.selectImagingTool(img=ks_mat_);
            soln.fileprefix = strrep(this.artery.fileprefix, "_radialarterylee", "_radialarteryleeks");
            this.product_ = soln;      

            % plot idealized artery
            h = solved_star.plot( ...
                tag=this.artery.fileprefix, ...
                xlim=[-10, this.times_sampled(end)+10], ...
                zoomMeas=1, ...
                zoomModel=1);
            saveFigure2(h, ...
                this.artery.fqfp + "_" + stackstr(), ...
                closeFigure=this.closeFigures);
        end
        function [qs,A_qs,qs_] = build_simulated(this, ks, times_sampled)
            %% BUILD_SIMULATED simulates tissue activity with passed and internal parameters.
            %  ks double is [k1 k2 k3 k4 k5 k6 k7 k8 k9] ~ 
            %               {'\alpha' '\beta' 'p' 'dp_2' 't_0' 'steadystate\_fraction' 'recirc\_fraction' 'recirc\_delay' 'baseline\_fraction'}.

            arguments
                this mlaif.ArteryLee2021Model
                ks {mustBeNumeric}
                times_sampled {mustBeNumeric} = []
            end

            [qs,A_qs,qs_]  = this.sampled(ks, this.Data, [], times_sampled);
        end
        function rho = deconvolved(this)
            rho = this.solver_.ArteryInterpolated;
        end
        function Q = loss(this)
            Q = this.solver_.loss();
        end
        function h = plot(this)
            N = length(this.kernel);
            z = this.solver_.M0*trapz(this.solver_.ArteryInterpolated)/trapz(this.kernel);
            h = this.solver_.plot( ...
                activityUnits="Bq", ...
                xlim=[-10, 200], ...
                xs={0:N-1}, ys={this.kernel}, zooms={z});
        end
        function writetable(this, fqfileprefix, opts)
            arguments
                this mlswisstrace.RadialArteryLee2021
                fqfileprefix {mustBeTextScalar} = this.fqfileprefix
                opts.scaling double = 1
            end

            this.fqfileprefix = fqfileprefix;
            rho_Bq_mL = opts.scaling*this.deconvolved()';
            N = length(rho_Bq_mL);
            times_sec = (0:N-1)';
            T = table(times_sec, rho_Bq_mL);
            writetable(T, this.fqfileprefix+stackstr(3)+".csv");
        end
    end
    
    methods (Static)
        function this = create(opts)
            arguments
                opts.artery mlfourd.ImagingContext2
                opts.closeFigures logical = false
                opts.kernel {mustBeNumeric}
                opts.model_kind {mustBeTextScalar} = "3bolus"
                opts.Nensemble double = 10
                opts.t0_forced {mustBeNumeric} = []
                opts.tracer {mustBeTextScalar} = "Unknown"
                opts.scalingWorldline {mustBeNumeric} = 1
            end

            this = mlaif.RadialArteryLee2024Model();

            this.artery_= mlpipeline.ImagingMediator.ensureFiniteImagingContext(opts.artery);
            ifc = this.artery.imagingFormat;
            this.measurement_ = asrow(ifc.img);
            this.times_sampled_ = asrow(ifc.json_metadata.timesMid);

            this.closeFigures = opts.closeFigures;
            this.Nensemble = opts.Nensemble;

            this.Data.kernel = opts.kernel;
            this.Data.model_kind = opts.model_kind;
            this.Data.N = floor(this.times_sampled(end)) + 1;
            if ~isempty(opts.t0_forced)
                this.Data.t0_forced = opts.t0_forced;
            end
            this.Data.tracer = opts.tracer;
            this.scalingWorldline = opts.scalingWorldline;

            this.LENK = 12;
        end

        function vec = apply_dispersion(vec_, Data)
            %% Args:            
            %      vec double, has 1 Hz sampling
            %      Data struct
            %  Returns:
            %      vec double, sampled at opts.Data.timesMid

            arguments
                vec_ double 
                Data struct = []
            end      
            if isempty(Data)
                return
            end
            vec = conv(vec_, Data.kernel);
            vec = vec(1:length(vec_));
        end
        function loss = loss_function(ks, kernel, tracer, model_kind, measurement_)
            import mlaif.RadialArteryLee2024Model.sampled
            import mlaif.ArteryLee2021Model.decay_corrected
            
            Data.N = length(measurement_);
            Data.tracer = tracer;
            Data.model_kind = model_kind;
            Data.kernel = kernel;

            estimation = sampled(ks, Data, [], []); % \in [0 1] 
            estimation = estimation/max(estimation);
            measurement_norm = measurement_/max(measurement_); % \in [0 1] 
            positive = measurement_norm > 0.01;
            eoverm = estimation(positive)./measurement_norm(positive);
            
            estimation_dc = decay_corrected(estimation, tracer, 0);
            estimation_dc_norm = estimation_dc/max(estimation_dc); % \in [0 1]  
            measurement_dc = decay_corrected(measurement_, tracer, 0);
            measurement_dc_norm = measurement_dc/max(measurement_dc); % \in [0 1]             
            eoverm_dc = estimation_dc_norm(positive)./measurement_dc_norm(positive);
            
            loss = mean(abs(1 - 0.5*eoverm - 0.5*eoverm_dc));
        end
        function m = preferredMap()
            m = mlaif.ArteryLee2021Model.preferredMap;
            m('k10') = struct('min',  0,     'max',  15,    'init', 10,    'sigma', 0.05); % recirc delay 
            m('k11') = struct('min',  0.5,   'max',   1,    'init',  0.95, 'sigma', 0.05); % amplitude fraction above baseline \approx 0.95
        end    
        function [qs,A_qs,qs_] = sampled(ks, Data, artery_interpolated, times_sampled)
            %% Returns:
            %      qs, the Bayesian estimate of the measured AIF, including baseline.
            %      qs_, the Bayesian estimate of the idealized AIF, scaled to unity.
            %      A_qs, providing the scaling factor lost to dispersion;
            
            arguments
                ks {mustBeNumeric}
                Data struct
                artery_interpolated {mustBeNumeric} = []
                times_sampled {mustBeNumeric} = []
            end

            amplitude = ks(11); % signal/(signal + baseline)
            baseline_frac = 1 - amplitude;
            qs_ = mlaif.RadialArteryLee2024Model.solution(ks, Data);             
            qs__ = mlaif.RadialArteryLee2024Model.apply_dispersion(qs_, Data); % 1 Hz sampling -> sampling of times_sampled
            A_qs = 1/max(qs__); % amplitude lost to dispersion > 1
            qs = amplitude*qs__ + baseline_frac*max(qs__);

            if ~isempty(times_sampled)
                idx_sampled = round(times_sampled - times_sampled(1)) + 1;
                qs = qs(idx_sampled);
                qs_ = qs_(idx_sampled);
            end                
        end
        function qs = solution(ks, Data)
            qs = mlaif.ArteryLee2021Model.solution(ks, Data);
        end
    end

    %% PROTECTED

    properties (Access = protected)
        artery_
        artery_interpolated_
        measurement_
        product_
        solver_ % for solver_solve
        times_sampled_
    end
    
    methods (Access = protected)
        function adjustMapForTracer(this)
            this.map = mlaif.ArteryLee2021Model.static_adjustMapForTracer(this.map, this.Data);
        end
        function adjustMapForModelKind(this)
            this.map = mlaif.ArteryLee2021Model.static_adjustMapForModelKind(this.map, this.Data);
        end
        function that = copyElement(this)
            %%  See also web(fullfile(docroot, 'matlab/ref/matlab.mixin.copyable-class.html'))
            
            that = copyElement@matlab.mixin.Copyable(this);
            that.solver_ = copy(this.solver_);
        end
        function img_int = interp1_artery(this, opts)
            arguments
                this mlaif.RadialArteryLee2024Model
                opts.times_sampled {mustBeNumeric} = this.times_sampled
                opts.img {mustBeNumeric} = this.artery.imagingFormat.img
                opts.timesAll {mustBeNumeric} = 0:this.times_sampled(end)
            end

            img_int = interp1( ...
                [0, opts.times_sampled], [0, opts.img], ...
                opts.timesAll, "linear", "extrap");
        end
        function solved = solver_solve(this)
            if isa(this.solver_, "mloptimization.SimulatedAnnealing")
                solved = this.solver_.solve(@mlaif.RadialArteryLee2024Model.loss_function);
                return
            end
            if isa(this.solver_, "mlnest.MultiNest")
                solved = this.solver_.solve(@mlaif.RadialArteryLee2024Model.signalmodel);
                return
            end
            error("mlaif:TypeError", stackstr())
        end

        function this = RadialArteryLee2024Model()
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
