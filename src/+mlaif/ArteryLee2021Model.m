classdef ArteryLee2021Model < handle & mlsystem.IHandle
	%% ArteryLee2021Model  

	%  $Revision$
 	%  was created 14-Mar-2021 17:21:03 by jjlee,
 	%  last modified $LastChangedDate$ and placed into repository /Users/jjlee/MATLAB-Drive/mlaif/src/+mlaif.
 	%% It was developed on Matlab 9.9.0.1592791 (R2020b) Update 5 for MACI64.  Copyright 2021 John Joowon Lee.
 	
    properties (Constant)
        ks_names = {'\alpha' '\beta' 'p' 'dp_F' 't_0' 'steadystate\_fraction' 'bolus_2\_fraction' 'bolus_2\_delay' 'recirc\_fraction' 'recirc\_delay' 'amplitude' '\gamma'}
    end
    
	properties 	
        closeFigures
        Data
        LENK
 		map  
        Nensemble
    end

    properties (Dependent)
        artery % mlfourd.ImagingContext2, entrypoint of arterial data
        artery_interpolated % double
        measurement % double
        model % mlaif.ArteryLee2021Model
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

            this.artery_interpolated_ = this.interp1_artery();
            g = this.artery_interpolated_;
        end
        function g = get.measurement(this)
            g = this.measurement_;
        end
        function g = get.model(this)
            g = this; % no copy for performance
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
                this mlaif.ArteryLee2021Model
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
                this.solver_ = mlpet.ArterySimulAnneal(context=this);
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
                solved{idx} = this.solver_.solve(@mlaif.ArteryLee2021Model.loss_function);
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
            img_new_ = solved_star.rescaleModelEstimate(A*ideal);
            ai_new_ = copy(this.artery.selectImagingTool(img=img_new_));
            ai_new_.fileprefix = strrep(this.artery.fileprefix, "_pet", "_arterylee");
            this.artery_ = ai_new_;  
            solved_star.ArteryInterpolated = this.interp1_artery(img=img_new_);   
            this.solver_ = solved_star; 

            % product_ := ks
            ks_mat_= [asrow(pr_.ks), pr_.loss];
            ks_mat_ = single(ks_mat_);
            soln = this.artery.selectImagingTool(img=ks_mat_);
            soln.fileprefix = strrep(this.artery.fileprefix, "_pet", "_arteryleeks");
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
        function [k,sk] = k1(this, varargin)
            [k,sk] = k1(this.solver_, varargin{:});
        end
        function [k,sk] = k2(this, varargin)
            [k,sk] = k2(this.solver_, varargin{:});
        end
        function [k,sk] = k3(this, varargin)
            [k,sk] = k3(this.solver_, varargin{:});
        end
        function [k,sk] = k4(this, varargin)
            [k,sk] = k4(this.solver_, varargin{:});
        end
        function [k,sk] = k5(this, varargin)
            [k,sk] = k5(this.solver_, varargin{:});
        end
        function [k,sk] = k6(this, varargin)
            [k,sk] = k6(this.solver_, varargin{:});
        end
        function [k,sk] = k7(this, varargin)
            [k,sk] = k7(this.solver_, varargin{:});
        end
        function [k,sk] = k8(this, varargin)
            [k,sk] = k8(this.solver_, varargin{:});
        end
        function [k,sk] = k9(this, varargin)
            [k,sk] = k9(this.solver_, varargin{:});
        end
        function [k,sk] = k10(this, varargin)
            [k,sk] = k10(this.solver_, varargin{:});
        end
        function [k,sk] = k11(this, varargin)
            [k,sk] = k11(this.solver_, varargin{:});
        end
        function [k,sk] = k12(this, varargin)
            [k,sk] = k12(this.solver_, varargin{:});
        end
        function [k,sk] = ks(this, varargin)
            k = zeros(1,this.LENK);
            sk = zeros(1,this.LENK);
            [k(1),sk(1)] = k1(this.solver_, varargin{:});
            [k(2),sk(2)] = k2(this.solver_, varargin{:});
            [k(3),sk(3)] = k3(this.solver_, varargin{:});
            [k(4),sk(4)] = k4(this.solver_, varargin{:});
            [k(5),sk(5)] = k5(this.solver_, varargin{:});
            [k(6),sk(6)] = k6(this.solver_, varargin{:});
            [k(7),sk(7)] = k7(this.solver_, varargin{:});
            [k(8),sk(8)] = k8(this.solver_, varargin{:});
            [k(9),sk(9)] = k9(this.solver_, varargin{:});
            [k(10),sk(10)] = k10(this.solver_, varargin{:});
            [k(11),sk(11)] = k11(this.solver_, varargin{:});
            [k(12),sk(12)] = k12(this.solver_, varargin{:});
        end	
        
        %% UTILITIES

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
    end 
    
    methods (Static)
        function this = create(opts)
            arguments
                opts.artery mlfourd.ImagingContext2
                opts.closeFigures logical = false
                opts.model_kind {mustBeTextScalar} = "3bolus"
                opts.Nensemble double = 1
                opts.t0_forced {mustBeNumeric} = []
                opts.tracer {mustBeTextScalar} = "Unknown"
            end

            this = mlaif.ArteryLee2021Model();

            this.artery_= mlpipeline.ImagingMediator.ensureFiniteImagingContext(opts.artery);
            ifc = this.artery.imagingFormat;
            this.measurement_ = asrow(ifc.img);
            this.times_sampled_ = asrow(ifc.json_metadata.timesMid);

            this.closeFigures = opts.closeFigures;
            this.Nensemble = opts.Nensemble;

            this.Data.model_kind = opts.model_kind;
            this.Data.N = floor(this.times_sampled(end)) + 1;
            if ~isempty(opts.t0_forced)
                this.Data.t0_forced = opts.t0_forced;
            end
            this.Data.tracer = opts.tracer;

            this.LENK = 12;
        end

        function loss = loss_function(ks, Data, ~, times_sampled, measurement)
            import mlaif.ArteryLee2021Model.sampled
            import mlaif.ArteryLee2021Model.decay_corrected
            
            estimation = sampled(ks, Data, [], times_sampled); % \in [0 1] 
            estimation = estimation/max(estimation);
            measurement_norm = measurement/max(measurement); % \in [0 1] 
            positive = measurement_norm > 0.01;
            eoverm = estimation(positive)./measurement_norm(positive);
            
            estimation_dc = decay_corrected(estimation, Data.tracer, 0);
            estimation_dc_norm = estimation_dc/max(estimation_dc); % \in [0 1]  
            measurement_dc = decay_corrected(measurement, Data, 0);
            measurement_dc_norm = measurement_dc/max(measurement_dc); % \in [0 1]             
            eoverm_dc = estimation_dc_norm(positive)./measurement_dc_norm(positive);
            
            loss = mean(abs(1 - 0.5*eoverm - 0.5*eoverm_dc));
        end
        function [qs,A_qs,qs_] = sampled(ks, Data, artery_interpolated, times_sampled)
            %% Returns:
            %      qs, the Bayesian estimate of the measured boxcar AIF, including baseline.
            %      qs_, the Bayesian estimate of the idealized AIF, scaled to unity.
            %      A_qs, providing the scaling factor lost to the boxcar.
            
            arguments
                ks {mustBeNumeric}
                Data struct
                artery_interpolated {mustBeNumeric} = []
                times_sampled {mustBeNumeric} = []
            end
            
            amplitude = ks(11);
            qs_ = amplitude*mlaif.ArteryLee2021Model.solution(ks, Data); % \in [0 1]
            qs = qs_;
            A_qs = 1/max(qs);

            if ~isempty(times_sampled)
                idx_sampled = round(times_sampled - times_sampled(1)) + 1;
                qs = qs(idx_sampled);
                qs_ = qs_(idx_sampled);
            end
        end
        function qs = solution(ks, Data)
            %% @return the idealized true AIF without baseline, scaled to unity.
            
            N = Data.N;
            tracer = Data.tracer;
            model_kind = Data.model_kind;
            import mlaif.ArteryLee2021Model
            switch model_kind
                case '1bolus'
                    qs = ArteryLee2021Model.solution_1bolus(ks, N, tracer, ks(3));
                case '2bolus'
                    qs = ArteryLee2021Model.solution_2bolus(ks, N, tracer, ks(3));
                case '3bolus'
                    qs = ArteryLee2021Model.solution_3bolus(ks, N, tracer);
                case '4bolus'
                    qs = ArteryLee2021Model.solution_4bolus(ks, N, tracer);
                otherwise
                    error('mlaif:ValueError', ...
                        'ArteryLee2021Model.solution.model_kind = %s', model_kind)
            end
        end
        function qs = solution_1bolus(ks, N, tracer, p)
            %% stretched gamma distribution

            import mlaif.ArteryLee2021Model.slide
            import mlaif.ArteryLee2021Model.halflife
            t = 0:N-1;
            t0 = ks(5);
            a = ks(1);
            b = ks(2);
            
            if (t(1) >= t0) 
                t_ = t - t0;
                qs = t_.^a .* exp(-(b*t_).^p);
            else % k is complex for t - t0 < 0
                t_ = t - t(1);
                qs = t_.^a .* exp(-(b*t_).^p);
                qs = slide(qs, t, t0 - t(1));
            end
            assert(all(imag(qs) == 0))
            % if ~isempty(tracer)
            %     qs = qs .* 2.^(-t/halflife(tracer));
            % end
            qs = qs/max(qs); % \in [0 1] 
        end
        function qs = solution_2bolus(ks, N, tracer, p)
            %% stretched gamma distribution + steadystate

            import mlaif.ArteryLee2021Model.slide
            import mlaif.ArteryLee2021Model.halflife
            t = 0:N-1;
            t0 = ks(5);
            a = ks(1);
            b = ks(2);
            g = ks(12);
            ss_frac = ks(6);
            dp = ks(4);
            p2 = max(p + dp, 0);
            
            if (t(1) >= t0) 
                t_ = t - t0;
                k_ = t_.^a .* exp(-(b*t_).^p);
                k_ = k_/max(k_); % else overwhelmed by ss_
                [~,ttp_] = max(k_);
                ss_ = (1 - exp(-t_/ttp_)).*exp(-(g*t_).^p2);
                ss_ = ss_/max(ss_);
                qs = (1 - ss_frac)*k_ + ss_frac*ss_;
            else % k is complex for t - t0 < 0
                t_ = t - t(1);
                k_ = t_.^a .* exp(-(b*t_).^p);
                k_ = k_/max(k_); % else overwhelmed by ss_
                [~,ttp_] = max(k_);
                ss_ = (1 - exp(-t_/ttp_)).*exp(-(g*t_).^p2);
                qs = (1 - ss_frac)*k_ + ss_frac*ss_;
                qs = slide(qs, t, t0 - t(1));
            end
            assert(all(imag(qs) == 0))
            % if ~isempty(tracer)
            %     qs = qs .* 2.^(-t/halflife(tracer));
            % end
            qs = qs/max(qs); % \in [0 1] 
        end
        function qs = solution_3bolus(ks, N, tracer)
            %% stretched gamma distribution + recirc stretched gamma distribution + rising steadystate; 
            %  forcing p2 = p - dp2 < p, to be more dispersive

            import mlaif.ArteryLee2021Model.solution_1bolus
            import mlaif.ArteryLee2021Model.solution_2bolus
            import mlaif.ArteryLee2021Model.slide
            recirc_frac = ks(9);
            recirc_delay = ks(10);
            
            qs1 = solution_2bolus(ks, N, tracer, ks(3));
            qs_recirc = solution_1bolus(ks, N, tracer, ks(3) + ks(4));
            qs_recirc = slide(qs_recirc, 0:N-1, recirc_delay);
            qs = (1 - recirc_frac)*qs1 + recirc_frac*qs_recirc;
            qs = qs/max(qs); % \in [0 1] 
        end
        function qs = solution_4bolus(ks, N, tracer)
            %% 2x stretched gamma distribution + recirc stretched gamma distribution + rising steadystate; 
            %  forcing p2 = p - dp2 < p, to be more dispersive

            import mlaif.ArteryLee2021Model.solution_3bolus
            import mlaif.ArteryLee2021Model.slide
            bolus2_frac = ks(7);
            bolus2_delay = ks(8);
            
            qs1 = solution_3bolus(ks, N, tracer);
            qs1 = qs1/max(qs1);
            qs2 = solution_3bolus(ks, N, tracer);
            qs2 = slide(qs2, 0:N-1, bolus2_delay);
            qs2 = qs2/max(qs2);     

            qs = (1 - bolus2_frac)*qs1 + bolus2_frac*qs2;
            qs = qs/max(qs); % \in [0 1] 
        end
        
        %% UTILITIES

        function rho = decay_corrected(rho, tracer, t0)
            arguments
                rho double
                tracer {mustBeTextScalar}
                t0 double = 0
            end
            import mlaif.ArteryLee2021Model.halflife
            times = (0:(length(rho)-1)) - t0;
            times(1:t0+1) = 0;
            rho = rho .* 2.^(times/halflife(tracer));
        end 
        function [vec,T] = ensureRow(vec)
            if (~isrow(vec))
                vec = vec';
                T = true;
                return
            end
            T = false; 
        end
        function tau = halflife(tracer)
            switch convertStringsToChars(upper(tracer))
                case {'FDG' '18F'}
                    tau = 1.82951 * 3600; % +/- 0.00034 h * sec/h
                case {'RO948', 'MK6240', 'GTP1', 'ASEM', 'AZAN'}
                    tau = 1.82951 * 3600; % +/- 0.00034 h * sec/h
                case {'HO' 'CO' 'OC' 'OO' '15O' 'WATER' 'CARBONMONOXIDE' 'CARBON_MONOXIDE' 'OXYGEN'}
                    tau = 122.2416;
                otherwise
                    error('mlaif:ValueError', ...
                        'ArteryLee2021Model.halflife.tracer = %s', tracer)
            end            
        end
        function conc = slide(conc, t, Dt)
            %% SLIDE slides discretized function conc(t) to conc(t - Dt);
            %  Dt > 0 will slide conc(t) towards later times t.
            %  Dt < 0 will slide conc(t) towards earlier times t.
            %  It works for inhomogeneous t according to the ability of interp1 to interpolate.
            %  It may not preserve information according to the Nyquist-Shannon theorem.  
            
            import mlaif.ArteryLee2021Model;
            [conc,trans] = ArteryLee2021Model.ensureRow(conc);
            t            = ArteryLee2021Model.ensureRow(t);
            
            tspan = t(end) - t(1);
            tinc  = t(2) - t(1);
            t_    = [(t - tspan - tinc) t];   % prepend times
            conc_ = [zeros(size(conc)) conc]; % prepend zeros
            conc_(isnan(conc_)) = 0;
            conc  = interp1(t_, conc_, t - Dt); % interpolate onto t shifted by Dt; Dt > 0 shifts to right
            
            if (trans)
                conc = conc';
            end
        end

        %% maps

        function m = preferredMap()
            T = 120;  % sec
            m = containers.Map;
            m('k1')  = struct('min',  0.5,   'max',   2,    'init',  0.5,  'sigma', 0.05); % alpha
            m('k2')  = struct('min',  0.05,  'max',   0.15, 'init',  0.05, 'sigma', 0.05); % beta in 1/sec
            m('k3')  = struct('min',  0.25,  'max',   3,    'init',  1,    'sigma', 0.05); % p
            m('k4')  = struct('min', -1,     'max',   0,    'init', -0.5,  'sigma', 0.05); % dp steady-state; p_ss = max(p + dp, 0)
            m('k5')  = struct('min',  5,     'max',  30,    'init',  0,    'sigma', 0.05); % t0 in sec
            m('k6')  = struct('min',  0,     'max',   0.2,  'init',  0.05, 'sigma', 0.05); % steady-state fraction in (0, 1), for rising baseline
            m('k7')  = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % 2nd bolus fraction
            m('k8')  = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % 2nd bolus delay
            m('k9')  = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % recirc fraction 
            m('k10') = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % recirc delay 
            m('k11') = struct('min',  1,     'max',   2,    'init',  1.3,  'sigma', 0.05); % amplitude 
            m('k12') = struct('min',  1/T,   'max',   T,    'init',  1/T,  'sigma', 0.05); % g for steady-state accumulation 
        end    

        function map = static_adjustMapForModelKind(map, Data)
            switch convertStringsToChars(lower(Data.model_kind))
                case '1bolus'        
                    map('k4') = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % dp steady-state; p_ss = max(p + dp, 0)

                    map('k6') = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % steady-state fraction in (0, 1), for rising baseline
                    map('k7') = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % 2nd bolus fraction
                    map('k8') = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % 2nd bolus delay
                    map('k9') = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % recirc fraction 
                    map('k10') = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % recirc delay 
                    map('k12') = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % g for steady-state accumulation 
                case '2bolus'
                    map('k4') = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % dp steady-state; p_ss = max(p + dp, 0)

                    map('k7') = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % 2nd bolus fraction
                    map('k8') = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % 2nd bolus delay
                    map('k9') = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % recirc fraction 
                    map('k10') = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % recirc delay 
                case '3bolus'
                    map('k7') = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % 2nd bolus fraction
                    map('k8') = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % 2nd bolus delay
                case '4bolus'
                otherwise
            end

            if isfield(Data, "t0_forced") && ~isempty(Data.t0_forced)
                t0_forced = Data.t0_forced;
                map('k5') = struct('min', t0_forced-1, 'max', t0_forced+1, 'init', t0_forced, 'sigma', 0.05); % t0 in sec
            end
        end

        function map = static_adjustMapForTracer(map, Data)
            switch convertStringsToChars(upper(Data.tracer))
                case {'RO948', 'MK6240', 'GTP1', 'ASEM', 'AZAN', 'VAT'}
                    T = 7200;  % sec
                    map('k6')  = struct('min',  0.05,  'max',   0.5,  'init',  0.05, 'sigma', 0.05); % steady-state fraction in (0, 1)  
                    map('k9')  = struct('min',  0.05,  'max',   0.25, 'init',  0.05, 'sigma', 0.05); % recirc fraction < 0.5, for 2nd bolus
                    map('k10') = struct('min',  5,     'max',  20,    'init', 10,    'sigma', 0.05); % recirc delay in sec
                    map('k12') = struct('min',  1/T,   'max', 500/T,  'init',  1/T,  'sigma', 0.05); % g for steady-state accumulation
                    % map('k1')  = struct('min',  1,     'max',   8,    'init',  5,    'sigma', 0.05); % alpha
                    % map('k2')  = struct('min',  1e-2,  'max',   0.1,  'init',  0.1,  'sigma', 0.05); % beta
                    % map('k3')  = struct('min',  0.5,   'max',   1,    'init',  0.5,  'sigma', 0.05); % p
                    % map('k4')  = struct('min', -0.5,   'max',   0,    'init',  0,    'sigma', 0.05); % dp steady-state; p_ss = max(p + dp, 0)
                    % map('k5')  = struct('min',  5,     'max',  70,    'init',  0,    'sigma', 0.05); % t0 in sec
                    % map('k6')  = struct('min',  0,     'max',   0.5,  'init',  0,    'sigma', 0.05); % steady-state fraction in (0, 1)  
                    % map('k7')  = struct('min',  0,     'max',   0.2,  'init',  0,    'sigma', 0.05); % 2nd bolus fraction
                    % map('k8')  = struct('min',  0,     'max',  15,    'init',  0,    'sigma', 0.05); % 2nd bolus delay
                    % map('k9')  = struct('min',  0,     'max',   0.2,  'init',  0,    'sigma', 0.05); % recirc fraction < 0.5, for 2nd bolus
                    % map('k10') = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % recirc delay in sec
                    % map('k12') = struct('min',  0,     'max',   5e-3, 'init',  0,    'sigma', 0.05); % g for steady-state accumulation & 2nd bolus
                case {'FDG' '18F'} 
                    T = 3600;  % sec
                    map('k6')  = struct('min',  0.05,  'max',   0.5,  'init',  0.05, 'sigma', 0.05); % steady-state fraction in (0, 1)  
                    map('k9')  = struct('min',  0.05,  'max',   0.25, 'init',  0.05, 'sigma', 0.05); % recirc fraction < 0.5, for 2nd bolus
                    map('k10') = struct('min',  5,     'max',  20,    'init', 10,    'sigma', 0.05); % recirc delay in sec
                    map('k12') = struct('min',  1/T,   'max', 500/T,  'init',  1/T,  'sigma', 0.05); % g for steady-state accumulation
                case {'HO' 'OH'}
                case {'CO' 'OC' 'CARBONMONOXIDE' 'CARBON_MONOXIDE'}  
                    map('k1')  = struct('min',  0.1,  'max',  10,    'init',  0.25, 'sigma', 0.05); % alpha
                    map('k2')  = struct('min', 10,    'max',  25,    'init', 15,    'sigma', 0.05); % beta
                    map('k3')  = struct('min',  0.25, 'max',   1,    'init',  0.5,  'sigma', 0.05); % p
         
                    %map('k6') = struct('min',  0.05, 'max',   0.5,  'init',  0.05, 'sigma', 0.05); % steady-state fraction in (0, 1)      
                    %map('k10') = struct('min', 30,    'max', 120,    'init', 60,    'sigma', 0.05); % recirc delay in sec
                case {'OO' 'OXYGEN'}
                    % allow double inhalation

                    map('k1')  = struct('min',  2,     'max',  10,    'init',  0.05, 'sigma', 0.05); % alpha
                    map('k2')  = struct('min',  0.8,   'max',   2,    'init',  0.05, 'sigma', 0.05); % beta
                    map('k3')  = struct('min',  1,     'max',   3,    'init',  0.5,  'sigma', 0.05); % p
                otherwise
                    % noninformative
            end
        end
    end
    
    %% PROTECTED
    
    properties (Access = protected)
        artery_
        artery_interpolated_
        measurement_
        product_
        solver_
        times_sampled_
    end

    methods (Access = protected)   
 		function this = ArteryLee2021Model()	
        end
        function that = copyElement(this)
            that = copyElement@matlab.mixin.Copyable(this);
            if ~isempty(this.artery_)
                that.artery_ = copy(this.artery_); end
            if ~isempty(this.product_)
                that.product_ = copy(this.product_); end
        end

        %% UTILITIES

        function adjustMapForTracer(this)
            this.map = mlaif.ArteryLee2021Model.static_adjustMapForTracer(this.map, this.Data);
        end
        function adjustMapForModelKind(this)
            this.map = mlaif.ArteryLee2021Model.static_adjustMapForModelKind(this.map, this.Data);
        end
        function img_int = interp1_artery(this, opts)
            arguments
                this mlaif.ArteryLee2021Model
                opts.times_sampled {mustBeNumeric} = this.times_sampled
                opts.img {mustBeNumeric} = this.artery.imagingFormat.img
                opts.timesAll {mustBeNumeric} = 0:this.times_sampled(end)
            end

            img_int = interp1( ...
                [0, opts.times_sampled], [0, opts.img], ...
                opts.timesAll, "linear", "extrap");
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

