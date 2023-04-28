classdef RadialArtery < handle & mlio.AbstractHandleIO & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %% RadialArtery provides a strategy design pattern for inferring cerebral AIFs
    %  from measurements sampling the radial artery and an arbitrary model kernel for delay and dispersion.
    %  
    %  Created 28-Apr-2023 14:34:29 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlaif/src/+mlaif.
    %  Developed on Matlab 9.14.0.2239454 (R2023a) Update 1 for MACI64.  Copyright 2023 John J. Lee.
    
	properties  
        measurement % expose for performance when used strategies for solve
        model       %
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
        function this = RadialArtery(varargin)
            %% RADIALARTERY 
            %  @param Measurement is a table:  [times, ]activityDensities in [sec, ]Bq/mL.
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
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'Measurement', [], @istable);
            addParameter(ip, 'solver', 'simulanneal', @ischar);
            parse(ip, varargin{:});
            ipr = ip.Results;
            this.measurement = ipr.Measurement;
            vns = this.measurement.Properties.VariableNames;
            assert(strcmp('times', vns{1}))
            assert(strcmp('activityDensities', vns{end}))

 			this.model = mlaif.RadialArteryModel(varargin{:});
                        
            switch lower(ipr.solver)
                case 'simulanneal'
                    this.strategy_ = mlaif.RadialArterySimulAnneal( ...
                        'context', this, varargin{:});
                otherwise
                    error('mlaif:NotImplementedError', ...
                        'RadialArtery.ipr.solver->%s', ipr.solver)
            end
        end

        function rho = deconvolved(this)
            M0 = max(this.measurement.activityDensities);
            N = floor(this.measurement.times(end)) + 1;
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
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
