classdef ArteryModel
    %% ARTERYMODEL supports a strategy design pattern for inferring cerebral AIFs
    %  from measurements of dynamic imaging and an arbitrary model kernel for delay and dispersion.
    %  Measurements should be decay-corrected as typically done by reconstruction methods.
    %  
    %  Created 28-Apr-2023 14:36:00 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlaif/src/+mlaif.
    %  Developed on Matlab 9.14.0.2239454 (R2023a) Update 1 for MACI64.  Copyright 2023 John J. Lee.
    
    properties (Constant)
        knames = {'\alpha' '\beta' 'p' 'dp_2' 't_0' 'steadystate\_fraction' 'recirc\_fraction' 'recirc\_delay' 'baseline\_fraction' 'scale'}
    end
    
	properties 	
        kernel
 		map
        model_kind
        t0_forced
        tracer
    end

	methods		  
 		function this = ArteryModel(opts)
            %%  Args:
            %   opts.tracer {mustBeTextScalar} = '15O'
            %   opts.model_kind = '3bolus'
            %   opts.kernel double = 1
            %   opts.map = []
            %   opts.t0_forced = []
            
            arguments
                opts.tracer {mustBeTextScalar} = '15O'
                opts.model_kind = '3bolus'
                opts.kernel double = 1
                opts.map = []
                opts.t0_forced = []
            end
            if isempty(opts.map)
                opts.map = mlaif.ArteryModel.preferredMap();
            end

            this.tracer = upper(opts.tracer);
            this.model_kind = opts.model_kind;
            this.map = opts.map;
            this.t0_forced = opts.t0_forced;
            this = this.adjustMapForTracer();
            this = this.adjustMapForModelKind();
            this.kernel = opts.kernel;		
 		end
    end 

    methods (Static)
        function rho = decay_corrected_from_table(T, tracer, opts)
            arguments
                T table
                tracer {mustBeTextScalar}
                opts.timeShift double = 0
                opts.sign = 1
            end
            timeShift = floor(opts.timeShift);
            times = T.timesMid - timeShift;
            f = mlaif.ArteryModel.decayCorrectionFactors( ...
                tracer=tracer, timeShift=timeShift, timeForDecayCorrection=0, times=times);
            rho = T.activityDensity .* f.^opts.sign;
            rho = asrow(rho);
        end
        function f = decayCorrectionFactors(opts)
            %% DECAYCORRECTIONFACTORS increases observed emissions, removing effects of radiodecay.
            %  See also:  https://niftypet.readthedocs.io/en/latest/tutorials/corrqnt.html
            %  Args:
            %      opts.tracer {mustBeTextScalar}
            %      opts.timeShift double = 0
            %      opts.timeForDecayCorrection double = 0
            %      opts.times double = [], as row | as col
            %  Returns:
            %      f, a vector with same shape at this.times.
            
            arguments
                opts.tracer {mustBeTextScalar}
                opts.timeShift double = 0
                opts.timeForDecayCorrection double = 0
                opts.times double = [] % as row | as col
            end

            import mlaif.ArteryModel.halflife
            
            lambda = log(2)/halflife(opts.tracer);
            times1 = opts.times - opts.timeForDecayCorrection - opts.timeShift;
            taus = (times1(2:end) - times1(1:end-1));
            taus(end+1) = taus(end);
            f = lambda*taus ./ (exp(-lambda*times1).*(1 - exp(-lambda*taus)));
        end
        function qs = deconvolved(M0, ks, N, kernel, tracer, model_kind)
            baseline_frac = ks(9);
            qs = mlaif.ArteryModel.solution(ks, N, tracer, model_kind);            
            if ~isscalar(kernel)
                % adjust kernel integral as trapz
                card_kernel = trapz(kernel);
                conv_sk = conv(qs, kernel); % ~1
                max_sk = max(conv_sk(1:N)); % ~1
                qs = qs*(card_kernel/max_sk); 
                qs = M0*((1 - baseline_frac)*qs + baseline_frac);
            end
            if contains(model_kind, "window")
                scale = ks(10);
                qs = scale*M0*((1 - baseline_frac)*qs + baseline_frac);
            end
        end
        function tau = halflife(tracer)
            arguments
                tracer char {mustBeTextScalar} = '18F'
            end

            switch upper(tracer)
                case {'FDG' '18F'}
                    tau = 1.82951 * 3600; % +/- 0.00034 h * sec/h
                case {'HO' 'CO' 'OC' 'OO' '15O'}
                    tau = 122.2416;
                otherwise
                    error('mlan:ValueError', ...
                        'ArteryModel.halflife.tracer = %s', tracer)
            end            
        end
        function p = jeffreys_prior(t)
            arguments
                t double
            end
            tmax = max(t);
            tmin = max(eps, min(t));
            p = 1./(t*log(tmax/tmin));
            p = p/sum(p);
        end
        function loss = loss_function(ks, kernel, tracer, model_kind, Measurement)
            %% generally needs matching to DC and DUC; consider also using Jeffreys' prior to weight times
            %  Jeffreys' prior:
            %    indices = 1:N;
            %    indices = indices(positive);
            %    p = mlaif.ArteryModel.jeffreys_prior(indices);

            import mlaif.ArteryModel.sampled
            import mlaif.ArteryModel.decay_corrected_from_table
            import mlaif.ArteryModel.jeffreys_prior

            measurement = asrow(Measurement.activityDensity); % discrete
            estimation = sampled(max(measurement), ks, length(measurement), kernel, tracer, model_kind);
            positive = measurement > 0.001*max(measurement);
            eoverm = estimation(positive)./measurement(positive);

            measurement_duc = decay_corrected_from_table(Measurement, tracer, timeShift=ks(5), sign=1);
            T_est = Measurement; T_est.activityDensity = ascol(estimation);
            estimation_duc = decay_corrected_from_table(T_est, tracer, timeShift=ks(5), sign=1);                    
            eoverm_duc = estimation_duc(positive)./measurement_duc(positive);            

            % timesMid = asrow(Measurement.timesMid);
            % timesMid = timesMid(positive);
            % eoverm_jeff = 0.5*jeffreys_prior(timesMid).*(eoverm + eoverm_duc);
            % eoverm_jeff = (sum(eoverm)/sum(eoverm_jeff))*eoverm_jeff;
            % loss = mean(abs(1 - eoverm_jeff));
            loss = mean(abs(1 - 0.5*(eoverm + eoverm_duc)));
            loss = double(loss);
        end
        function m = preferredMap()
            m = containers.Map;
            m('k1')  = struct('min',  0.5,   'max',   1.5,  'init',  0.5,  'sigma', 0.05); % alpha
            m('k2')  = struct('min',  0.05,  'max',   0.15, 'init',  0.05, 'sigma', 0.05); % beta in 1/sec
            m('k3')  = struct('min',  1,     'max',   3,    'init',  1,    'sigma', 0.05); % p
            m('k4')  = struct('min',  0,     'max',   1,    'init',  0,    'sigma', 0.05); % |dp2| for 2nd bolus
            m('k5')  = struct('min',  0,     'max',  30,    'init',  0,    'sigma', 0.05); % t0 in sec
            m('k6')  = struct('min',  0.1,   'max',   0.3,  'init',  0.2,  'sigma', 0.05); % steady-state fraction in (0, 1), for rising baseline
            m('k7')  = struct('min',  0,     'max',   0.1,  'init',  0.05, 'sigma', 0.05); % recirc fraction < 0.5, for 2nd bolus
            m('k8')  = struct('min', 15,     'max',  90,    'init', 30,    'sigma', 0.05); % recirc delay in sec
            m('k9')  = struct('min',  0,     'max',   0.02, 'init',  0,    'sigma', 0.05); % baseline amplitude fraction \approx 0.05
            m('k10') = struct('min',  1,     'max',   3,    'init',  1,    'sigma', 0.05); % scaling that is not conserved by sampled()
        end    
        function qs = sampled(M0, ks, N, kernel, tracer, model_kind)
            %% @return the Bayesian estimate of the measured AIF, including baseline, scaled to unity.
            %  no baseline modelled for model_kind ~ "window"
            
            baseline_frac = ks(9);            
            if ~isscalar(kernel)
                qs = mlaif.ArteryModel.solution(ks, N, tracer, model_kind);
                qs = conv(qs, kernel);
                qs = qs(1:N);
                qs = qs/max(qs);
                qs = M0*((1 - baseline_frac)*qs + baseline_frac); 
                return
            end     
            if contains(model_kind, "window")
                scale = ks(10);
                re = regexp(model_kind, "\S+window(?<W>\d+)\S*", "names");
                W = str2double(re.W);
                qs = mlaif.ArteryModel.solution(ks, N+W-1, tracer, model_kind);
                qs = scale*M0*((1 - baseline_frac)*qs + baseline_frac);
                qs = mlaif.ArteryModel.move_window(qs, W=W);
            end  
        end
        function qs = solution(ks, N, tracer, model_kind)
            %% @return the idealized true AIF without baseline, scaled to unity.
            
            import mlaif.ArteryModel
            if contains(model_kind, '1bolus')
                qs = ArteryModel.solution_1bolus(ks, N, tracer, ks(3));
            end
            if contains(model_kind, '2bolus')
                qs = ArteryModel.solution_2bolus(ks, N, tracer, ks(3));
            end
            if contains(model_kind, '3bolus')
                qs = ArteryModel.solution_3bolus(ks, N, tracer);
            end
        end
        function qs = solution_1bolus(ks, N, ~, p)
            %% stretched gamma distribution for decay-uncorrected

            import mlaif.ArteryModel.slide
            import mlaif.ArteryModel.halflife
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
            qs = qs .* 2.^(-t/halflife(tracer));
            qs = qs/max(qs); % \in [0 1] 
        end
        function qs = solution_2bolus(ks, N, ~, p)
            %% stretched gamma distribution + rising steadystate for decay-uncorrected

            import mlaif.ArteryModel.slide
            import mlaif.ArteryModel.halflife
            t = 0:N-1;
            t0 = ks(5);
            a = ks(1);
            b = ks(2);
            g = ks(2);
            ss_frac = ks(6);
            
            if (t(1) >= t0) 
                t_ = t - t0;
                k_ = t_.^a .* exp(-(b*t_).^p);
                ss_ = 1 - exp(-g*t_);
                qs = (1 - ss_frac)*k_ + ss_frac*ss_;
            else % k is complex for t - t0 < 0
                t_ = t - t(1);
                k_ = t_.^a .* exp(-(b*t_).^p);
                ss_ = 1 - exp(-g*t_);
                qs = (1 - ss_frac)*k_ + ss_frac*ss_;
                qs = slide(qs, t, t0 - t(1));
            end
            assert(all(imag(qs) == 0))
            qs = qs .* 2.^(-t/halflife(tracer));
            qs = qs/max(qs); % \in [0 1] 
        end
        function qs = solution_3bolus(ks, N, tracer)
            %% stretched gamma distribution + rising steadystate + auxiliary stretched gamma distribution 
            %  for decay-uncorrected;
            %  forcing p2 = p - dp2 < p, to be more dispersive

            import mlaif.ArteryModel.solution_1bolus
            import mlaif.ArteryModel.solution_2bolus
            import mlaif.ArteryModel.slide
            recirc_frac = ks(7);
            recirc_delay = ks(8);
            
            qs2 = solution_2bolus(ks, N, tracer, ks(3));
            qs1 = solution_1bolus(ks, N, tracer, ks(3) - ks(4));
            qs1 = slide(qs1, 0:N-1, recirc_delay);
            qs = (1 - recirc_frac)*qs2 + recirc_frac*qs1;
            qs = qs/max(qs); % \in [0 1] 
        end
        
        %% UTILITIES
          
        function [vec,T] = ensureRow(vec)
            if (~isrow(vec))
                vec = vec';
                T = true;
                return
            end
            T = false; 
        end
        function qs = move_window(qs, opts)
            arguments
                qs double % asrow in; as row out
                opts.W double % window length
            end
            
            P = length(qs)-opts.W+1;
            N = P/opts.W;
            P1 = P + opts.W - 1;
            L = tril(ones(P,P1), opts.W-1) - tril(ones(P,P1), -1);
            L = L/opts.W;
            qs = asrow(L*ascol(qs));
        end
        function conc = slide(conc, t, Dt)
            %% SLIDE slides discretized function conc(t) to conc(t - Dt);
            %  Dt > 0 will slide conc(t) towards later timesMid t.
            %  Dt < 0 will slide conc(t) towards earlier timesMid t.
            %  It works for inhomogeneous t according to the ability of interp1 to interpolate.
            %  It may not preserve information according to the Nyquist-Shannon theorem.  
            
            import mlaif.ArteryModel;
            [conc,trans] = asrow(conc);
            t            = asrow(t);
            
            tspan = t(end) - t(1);
            tinc  = t(2) - t(1);
            t_    = [(t - tspan - tinc) t];   % prepend timesMid
            conc_ = [zeros(size(conc)) conc]; % prepend zeros
            conc_(isnan(conc_)) = 0;
            conc  = interp1(t_, conc_, t - Dt); % interpolate onto t shifted by Dt; Dt > 0 shifts to right
            
            if (trans)
                conc = conc';
            end
        end
    end
        
    %% PROTECTED
    
    methods (Access = protected)     
        function this = adjustMapForModelKind(this)
            if contains(this.model_kind, '1bolus', IgnoreCase=true)     
                this.map('k4') = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % dp2 for 2nd bolus

                this.map('k6') = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % steady-state fraction in (0, 1), for rising baseline
                this.map('k7') = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % recirc fraction < 0.5, for 2nd bolus
                this.map('k8') = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % recirc delay in sec
            end
            if contains(this.model_kind, '2bolus', IgnoreCase=true)
                this.map('k4') = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % dp2 for 2nd bolus

                this.map('k7') = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % recirc fraction < 0.5, for 2nd bolus
                this.map('k8') = struct('min',  0,     'max',   0,    'init',  0,    'sigma', 0.05); % recirc delay in sec
            end
            if ~isempty(this.t0_forced)
                this.map('k5') = struct('min', this.t0_forced, 'max', this.t0_forced, 'init', this.t0_forced, 'sigma', 0.05); % t0 in sec
            end
        end
        function this = adjustMapForTracer(this)
            switch upper(this.tracer)
                case {'FDG' '18F'}
                    this.map('k5') = struct('min',  0,    'max',  30,    'init',  0,    'sigma', 0.05); % t0 in sec   
                    this.map('k6') = struct('min', 0.05,  'max',   0.5,  'init',  0.05, 'sigma', 0.05); % steady-state fraction in (0, 1)  
                    this.map('k7') = struct('min', 0.05,  'max',   0.25, 'init',  0.05, 'sigma', 0.05); % recirc fraction < 0.5, for 2nd bolus
                    this.map('k8') = struct('min', 5,     'max',  20,    'init', 10,    'sigma', 0.05); % recirc delay in sec
                case {'HO' 'OH'}
                case {'CO' 'OC'}  
                    this.map('k1') = struct('min',  0.1,  'max',  10,    'init',  0.25, 'sigma', 0.05); % alpha
                    this.map('k2') = struct('min', 10,    'max',  25,    'init', 15,    'sigma', 0.05); % beta
                    this.map('k3') = struct('min',  0.25, 'max',   1,    'init',  0.5,  'sigma', 0.05); % p

                    this.map('k5') = struct('min',  0,    'max',  30,    'init',  0,    'sigma', 0.05); % t0 in sec            
                    this.map('k6') = struct('min',  0.05, 'max',   0.5,  'init',  0.05, 'sigma', 0.05); % steady-state fraction in (0, 1)      

                    this.map('k8') = struct('min', 30,    'max', 120,    'init', 60,    'sigma', 0.05); % recirc delay in sec
                case 'OO'
                    % allow double inhalation
                    this.map('k1') = struct('min', 1,     'max',   5,    'init',  1,    'sigma', 0.05); % alpha
                    this.map('k2') = struct('min', 0.1,   'max',  25,    'init',  0.1,  'sigma', 0.05); % beta

                    this.map('k5') = struct('min', 5,     'max',  30,    'init',  0,    'sigma', 0.05); % t0 in sec   

                    this.map('k7') = struct('min', 0.05,  'max',   0.49, 'init',  0.05, 'sigma', 0.05); % recirc fraction < 0.5, for 2nd bolus
                    this.map('k8') = struct('min', eps,   'max',  30,    'init', 15,    'sigma', 0.05); % recirc delay in sec
                otherwise
                    % noninformative
            end
            if ~isempty(this.t0_forced)
                this.map('k5') = struct('min', this.t0_forced, 'max', this.t0_forced, 'init', this.t0_forced, 'sigma', 0.05); % t0 in sec
            end
        end
    end

    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
