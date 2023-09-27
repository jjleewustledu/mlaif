classdef RadialArtery < handle & mlio.AbstractHandleIO & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %% RadialArtery provides a strategy design pattern for inferring cerebral AIFs
    %  from measurements sampling the radial artery and an arbitrary model kernel for delay and dispersion.
    %  
    %  Created 28-Apr-2023 14:34:29 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlaif/src/+mlaif.
    %  Developed on Matlab 9.14.0.2239454 (R2023a) Update 1 for MACI64.  Copyright 2023 John J. Lee.
    
	properties  
        Measurement % expose for performance when used strategies for solve
        model       %
        units
    end

    properties (Dependent)
        varNameDecayCorrected
        strategy
    end

	methods %% GET
        function g = get.strategy(this)
            g = this.strategy_;
        end
        function g = get.varNameDecayCorrected(this)
            vars = this.Measurement.Properties.VariableNames;
            indicator_dc = contains(vars, 'dc', IgnoreCase=true) | contains(vars, 'decay', IgnoreCase=true);
            if any(indicator_dc)
                g = vars{indicator_dc};
                return
            end
            g = [];
        end
    end

    methods
        function this = RadialArtery(Measurement, units, opts)
            %% RADIALARTERY 
            %  @param Measurement is a table:  [times, ]activityDensity in [sec, ]Bq/mL.
            %  @param solver is in {'simulanneal'}.

            %  for mlaif.RadialArteryModel: 
            %  @param tracer, passed to model.
            %  @param model_kind is char, e.g. '1bolus', '2bolus', '3bolus', passed to model.
            %  @param map is a containers.Map.  Default := RadialArteryModel.preferredMap.
 			%  @param kernel is numeric, default := 1.
            %  @param t0_forced is scalar, default empty, passed to model.
            %
            %  for mlaif.RadialArterySimulAnneal:
            %  @param context is mlaif.RadialArtery.
            %  @param fileprefix.

            arguments
                Measurement table
                units struct
                opts.solver string {mustBeTextScalar} = "simulanneal"
                opts.tracer {mustBeTextScalar} = 'FDG'
                opts.model_kind {mustBeTextScalar} = '3bolus'
                opts.kernel double = 1
            end
            
            this.Measurement = Measurement;
            this.units = units;
 			this.model = mlaif.RadialArteryModel( ...
                tracer=opts.tracer, model_kind=opts.model_kind, kernel=opts.kernel);
            this.standardize_ctor_args();
                        
            switch opts.solver
                case "simulanneal"
                    this.strategy_ = mlaif.RadialArterySimulAnneal(context=this);
                otherwise
                    error("mlaif:NotImplementedError", ...
                        "RadialArtery.ipr.solver->%s", ipr.solver)
            end
        end

        function rho = deconvolved(this)
            M0 = max(this.Measurement.activityDensity);
            N = floor(this.Measurement.times(end)) + 1;
            ks = this.strategy_.ks;
            mdl = this.model;
            rho = M0*this.model.deconvolved(ks, N, mdl.kernel, mdl.tracer, mdl.model_kind);
        end
        function Q = loss(this)
            Q = this.strategy_.loss();
        end
        function h = plot(this, varargin)
            h = this.strategy_.plot(varargin{:});
        end
        function h = plot_dc(this, varargin)
            h = this.strategy_.plot_dc(varargin{:});
        end
        function this = solve(this, varargin)
            %% @param required loss_function is function_handle.
            
            this.strategy_ = solve(this.strategy_, @mlaif.RadialArteryModel.loss_function);
        end
        function writetable(this, fqfn)
            writetable(this.strategy_, fqfn)
        end
    end
    
    %% PROTECTED
    
    properties (Access = protected)
        strategy_ % for solve
    end
    
    methods (Access = protected)
        function that = copyElement(this)
            %%  See also web(fullfile(docroot, 'matlab/ref/matlab.mixin.copyable-class.html'))
            
            that = copyElement@matlab.mixin.Copyable(this);
            that.strategy_ = copy(this.strategy_);
        end
        function standardize_ctor_args(this)
            if 1 == length(this.Measurement.Properties.VariableNames)
                % assume activityDensity is sampled at 1 Hz
                times = 0:size(this.Measurement,1)-1;
                this.Measurement = addvars(this.Measurement, times, 1);
            end
            switch lower(char(this.units.times)) % convert times to seconds
                case {'s', 'sec', 'second', 'seconds'}
                case {'m', 'min', 'minute', 'minutes'}
                    this.Measurement.times = round(60*this.Measurement.times);
                case {'h', 'hr', 'hour', 'hours'}
                    this.Measurement.times = round(3600*this.Measurement.times);
                otherwise
                    error("mlaif:ValueError", ...
                        "RadialArtery.standardize_ctor_args.this.units.times->%s", this.units.times)
            end
            if ~isempty(this.varNameDecayCorrected)
                % convert external decay corrections to decay uncorreted
                factor = 2.^(-this.Measurement.times/this.model.halflife);
                this.Measurement.activityDensity = this.Measurement.(this.varNameDecayCorrected).*factor;
                this.units.activityDensity = this.units.(this.varNameDecayCorrected);
            end
            switch lower(char(this.units.activityDensity)) % convert activityDensity to Bq/mL
                case {'bqml', 'bq/ml'}
                case {'kbqml', 'kbq/ml'}
                    this.Measurement.activityDensity = 1e3*this.Measurement.activityDensity;
                case {'mbqml', 'mbq/ml'}
                    this.Measurement.activityDensity = 1e6*this.Measurement.activityDensity;
                case {'uciml', 'uci/ml'}
                    this.Measurement.activityDensity = 37e3*this.Measurement.activityDensity;
                otherwise
                    error("mlaif:ValueError", ...
                        "RadialArtery.standardize_ctor_args.this.units.activityDensity->%s", this.units.activityDensity)
            end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
