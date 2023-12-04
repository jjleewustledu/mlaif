classdef Artery < handle & mlsystem.IHandle
    %% ARTERY provides a strategy design pattern for inferring cerebral AIFs
    %  from measurements of dynamic imaging and an arbitrary model kernel for delay and dispersion.
    %  Measurements should be decay-corrected as typically done by reconstruction methods.
    %
    %  Is is DEPRECATED.  Prefer mlaif.ArteryLee2021Model, mlpet.ArterySimulAnneal.
    %  
    %  Created 28-Apr-2023 14:34:29 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlaif/src/+mlaif.
    %  Developed on Matlab 9.14.0.2239454 (R2023a) Update 1 for MACI64.  Copyright 2023 John J. Lee.
    
	properties  
        Measurement % expose for performance when using strategies for solve
        model       %
        units
    end

    properties (Dependent)
        strategy
    end

	methods %% GET
        function g = get.strategy(this)
            g = this.strategy_;
        end
    end

    methods
        function this = Artery(Measurement, units, opts)
            %% ARTERY 
            %  @param Measurement is a table constructed from dynamic imaging that is decay-corrected:  
            %                     [timesMid, ]activityDensity in [sec, ]Bq/mL.
            %  @param solver is in {'simulanneal'}.

            %  for mlaif.ArteryModel: 
            %  @param tracer, passed to model.
            %  @param model_kind is char, e.g. '1bolus', '2bolus', '3bolus', '3bolus-window10', etc., passed to model.
 			%  @param kernel is numeric, default := 1.
            %  @param t0_forced is scalar, default empty, passed to model.
            %
            %  for mlaif.ArterySimulAnneal:
            %  @param context is mlaif.Artery.

            arguments
                Measurement table
                units struct = struct("timesMid", "sec", "activityDensity", "Bq/mL")
                opts.solver string {mustBeTextScalar} = "simulanneal"
                opts.tracer {mustBeTextScalar} = '15O'
                opts.model_kind {mustBeTextScalar} = '3bolus'
                opts.kernel double = 1
            end
            
            this.Measurement = Measurement;
            this.units = units;
 			this.model = mlaif.ArteryModel( ...
                tracer=opts.tracer, model_kind=opts.model_kind, kernel=opts.kernel);
            this.standardize_ctor_args();
                        
            switch opts.solver
                case "simulanneal"
                    this.strategy_ = mlpet.ArterySimulAnneal(context=this);
                otherwise
                    error("mlaif:NotImplementedError", "Artery.ipr.solver->%s", ipr.solver)
            end
        end

        function rho = deconvolved(this)
            M0 = max(this.Measurement.activityDensity);
            N = size(this.Measurement, 1);
            ks = this.strategy_.ks;
            mdl = this.model;
            rho = this.model.deconvolved(M0, ks, N, mdl.kernel, mdl.tracer, mdl.model_kind);
        end
        function Q = loss(this)
            Q = this.strategy_.loss();
        end
        function h = plot(this, varargin)
            h = this.strategy_.plot(varargin{:});
        end
        function this = solve(this, varargin)
            %% @param required loss_function is function_handle.
            
            this.strategy_ = solve(this.strategy_, @mlaif.ArteryModel.loss_function);
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
                timesMid = 0.5:size(this.Measurement,1)-0.5;
                this.Measurement = addvars(this.Measurement, timesMid, 1);
            end
            switch lower(char(this.units.timesMid)) % convert timesMid to seconds
                case {'s', 'sec', 'second', 'seconds'}
                case {'m', 'min', 'minute', 'minutes'}
                    this.Measurement.timesMid = round(60*this.Measurement.timesMid);
                case {'h', 'hr', 'hour', 'hours'}
                    this.Measurement.timesMid = round(3600*this.Measurement.timesMid);
                otherwise
                    error("mlaif:ValueError", ...
                        "Artery.standardize_ctor_args.this.units.timesMid->%s", this.units.timesMid)
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
                        "Artery.standardize_ctor_args.this.units.activityDensity->%s", this.units.activityDensity)
            end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
