classdef RadialArterySimulAnneal < mloptimization.SimulatedAnnealing
    %% line1
    %  line2
    %  
    %  Created 28-Apr-2023 14:35:30 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlaif/src/+mlaif.
    %  Developed on Matlab 9.14.0.2239454 (R2023a) Update 1 for MACI64.  Copyright 2023 John J. Lee.
    
	properties
        ks0
        ks_lower
        ks_upper
        quiet = false
        visualize = false
        visualize_anneal = false
        zoom = 1
 	end
    
	properties (Dependent) 
        kernel  
        ks
        ks_names
        model_kind
        tracer
    end

	methods %% GET        
        function g = get.kernel(this)
            g = this.model.kernel;
        end
        function g = get.ks(this)
            g = this.results_.ks;
        end
        function g = get.ks_names(this)
            g = this.model.knames;
        end
        function g = get.model_kind(this)
            g = this.model.model_kind;
        end
        function g = get.tracer(this)
            g = this.model.tracer;
        end
    end

    methods
        function fprintfModel(this)
            fprintf('RadialArterySimulAnnel:\n');  
            fprintf('%s %s:\n', this.tracer, this.model_kind); 
            for ky = 1:length(this.ks)
                fprintf('\t%s = %g\n', this.ks_names{ky}, this.ks(ky));
            end 
            fprintf('\tloss = %g\n', this.loss())
            keys = this.map.keys;
            for ky = 1:length(this.ks)
                fprintf('\tmap(''%s'') => %s\n', this.ks_names{ky}, ...
                    join(struct2str(this.map(keys{ky}), orientation='horz')));
            end
        end
        function Q = loss(this)
            Q = this.results_.sse;
        end
        function h = plot_dc(this, varargin)
            ip = inputParser;
            addParameter(ip, 'showKernel', false, @islogical)
            addParameter(ip, 'xlim', [-10 200], @isnumeric)
            addParameter(ip, 'ylim', [], @isnumeric)
            addParameter(ip, 'zoom', 1, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            this.zoom = ipr.zoom;
            deconvolved = @mlaif.RadialArteryModel.deconvolved;
            sampled = @mlaif.RadialArteryModel.sampled;
            decay_corrected = @mlaif.RadialArteryModel.decay_corrected;
            decay_corrected_from_table = @mlaif.RadialArteryModel.decay_corrected_from_table;
            M_activityDensities = decay_corrected_from_table(this.Measurement, this.tracer);
            N = this.Measurement.times(end) + 1;
            M0 = max(M_activityDensities);
            
            h = figure;
            samp = M0*sampled(this.ks, N, this.kernel, this.tracer, this.model_kind);
            samp = decay_corrected(samp, this.tracer, 1:N);
            deconvolved = M0*deconvolved(this.ks, N, this.kernel, this.tracer, this.model_kind);
            times = 0:N-1;
            
            if isempty(this.zoom)
                this.zoom = max(deconvolved)/max(this.kernel)/2;
            end
            if this.zoom ~= 1
                leg_kern = sprintf('kernel x%g', this.zoom);
            else
                leg_kern = 'kernel';
            end
            if ipr.showKernel
                plot(asrow(this.Measurement.times), asrow(M_activityDensities), 'o', ...
                    times, samp, ':', ...
                    times, deconvolved, '-', ...
                    times, this.zoom*this.kernel, '--', 'LineWidth', 2)
                legend({'measured', 'estimated', 'deconvolved', leg_kern})
            else
                plot(asrow(this.Measurement.times), asrow(M_activityDensities), 'o', ...
                    times, samp, ':', 'LineWidth', 2)
                legend({'measured', 'estimated'}, 'FontSize', 10)
            end
            if ~isempty(ipr.xlim); xlim(ipr.xlim); end
            if ~isempty(ipr.ylim); ylim(ipr.ylim); end
            xlabel('times / s')
            ylabel('activity / (Bq/mL)')
            annotation('textbox', [.25 .5 .5 .2], 'String', sprintfModel(this), 'FitBoxToText', 'on', 'FontSize', 10, 'LineStyle', 'none')
            title([clientname(false, 2) ' DECAY-CORRECTED for ' this.tracer], 'Interpreter', 'none')
        end
        function h = plot(this, varargin)
            ip = inputParser;
            addParameter(ip, 'showKernel', false, @islogical)
            addParameter(ip, 'xlim', [-10 200], @isnumeric)
            addParameter(ip, 'ylim', [], @isnumeric)
            addParameter(ip, 'zoom', 1, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            this.zoom = ipr.zoom;
            M = this.Measurement;
            N = M.times(end) + 1;
            M0 = max(M.activityDensities);
                        
            h = figure;
            samp = M0*this.model.sampled(this.ks, N, this.kernel, this.tracer, this.model_kind);
            deconvolved = M0*this.model.deconvolved(this.ks, N, this.kernel, this.tracer, this.model_kind);
            times = 0:N-1;
            
            if isempty(this.zoom)
                this.zoom = max(deconvolved)/max(this.kernel)/2;
            end
            if this.zoom ~= 1
                leg_kern = sprintf('kernel x%g', this.zoom);
            else
                leg_kern = 'kernel';
            end
            if ipr.showKernel
                plot(asrow(M.times), asrow(M.activityDensities), 'o', ...
                    times, samp, ':', ...
                    times, deconvolved, '-', ...
                    times, this.zoom*this.kernel, '--', 'LineWidth', 2)
                legend({'measured', 'estimated', 'deconvolved', leg_kern})
            else
                plot(asrow(M.times), asrow(M.activityDensities), 'o', ...
                    times, samp, ':', 'LineWidth', 2)
                legend({'measured', 'estimated'}, 'FontSize', 10)
            end
            if ~isempty(ipr.xlim); xlim(ipr.xlim); end
            if ~isempty(ipr.ylim); ylim(ipr.ylim); end
            xlabel('times / s')
            ylabel('activity / (Bq/mL)')
            annotation('textbox', [.25 .5 .3 .3], 'String', sprintfModel(this), 'FitBoxToText', 'on', 'FontSize', 10, 'LineStyle', 'none')
            title(clientname(false, 2))
        end
        function save(this)
            save([this.fileprefix '.mat'], this);
        end
        function saveas(this, fn)
            save(fn, this);
        end
        function this = solve(this, varargin)
            %% @param required loss_function is function_handle.
            
            ip = inputParser;
            addRequired(ip, 'loss_function', @(x) isa(x, 'function_handle'))
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            options_fmincon = optimoptions('fmincon', ...
                'FunctionTolerance', 1e-12, ...
                'OptimalityTolerance', 1e-12, ...
                'TolCon', 1e-14, ...
                'TolX', 1e-14);
            if this.visualize_anneal
                options = optimoptions('simulannealbnd', ...
                    'AnnealingFcn', 'annealingboltz', ...
                    'FunctionTolerance', eps, ...
                    'HybridFcn', {@fmincon, options_fmincon}, ...
                    'InitialTemperature', 20, ...
                    'MaxFunEvals', 50000, ...
                    'ReannealInterval', 200, ...
                    'TemperatureFcn', 'temperatureexp', ...
                    'Display', 'diagnose', ...
                    'PlotFcns', {@saplotbestx,@saplotbestf,@saplotx,@saplotf,@saplotstopping,@saplottemperature});
            else
                options = optimoptions('simulannealbnd', ...
                    'AnnealingFcn', 'annealingboltz', ...
                    'FunctionTolerance', eps, ...
                    'HybridFcn', {@fmincon, options_fmincon}, ...
                    'InitialTemperature', 20, ...
                    'MaxFunEvals', 50000, ...
                    'ReannealInterval', 200, ...
                    'TemperatureFcn', 'temperatureexp');
            end
 			[ks_,sse,exitflag,output] = simulannealbnd( ...
                @(ks__) ipr.loss_function( ...
                       ks__, this.kernel, this.tracer, this.model_kind, this.Measurement), ...
                this.ks0, this.ks_lower, this.ks_upper, options); 
            
            this.results_ = struct('ks0', this.ks0, 'ks', ks_, 'sse', sse, 'exitflag', exitflag, 'output', output); 
            if ~this.quiet
                fprintfModel(this)
            end
            if this.visualize
                plot(this)
            end
        end 
        function s = sprintfModel(this)
            s = sprintf('RadialArterySimulAnnel:\n');
            s = [s sprintf('%s %s:\n', this.tracer, this.model_kind)];
            for ky = 1:length(this.ks)
                s = [s sprintf('\t%s = %g\n', this.ks_names{ky}, this.ks(ky))]; %#ok<AGROW>
            end
            s = [s sprintf('\tloss = %g\n', this.loss())];
            keys = this.map.keys;
            for ky = 1:length(this.ks)
                s = [s sprintf('\tmap(''%s'') => %s\n', this.ks_names{ky}, ...
                    join(struct2str(this.map(keys{ky}), orientation='horz')))]; %#ok<AGROW>
            end
        end
        function writetable(this, fqfn)
            M0 = max(this.Measurement.activityDensities);
            N = floor(this.Measurement.times(end)) + 1;
            activityDensitiesDC = ascol(M0*this.model.sampled(this.ks, N, this.kernel, this.tracer, this.model_kind));
            times = (0:length(activityDensitiesDC)-1)'; % sec
            activityDensitiesDC = activityDensitiesDC .* 2.^(times/6582);
            activityDensitiesDC = activityDensitiesDC/37000;
            T = table(times, activityDensitiesDC);
            writetable(T, fqfn);
        end

 		function this = RadialArterySimulAnneal(varargin)
            this = this@mloptimization.SimulatedAnnealing(varargin{:}); 			
            [this.ks_lower,this.ks_upper,this.ks0] = remapper(this);
 		end
    end
    
    %% PROTECTED
    
    methods (Access = protected)
        function [m,sd] = find_result(this, lbl)
            ks_ = this.ks;
            assert(strcmp(lbl(1), 'k'))
            ik = str2double(lbl(2));
            m = ks_(ik);
            sd = 0;
        end
        function [lb,ub,ks0] = remapper(this)
            for i = 1:this.map.Count
                lbl = sprintf('k%i', i);
                lb(i)  = this.map(lbl).min; %#ok<AGROW>
                ub(i)  = this.map(lbl).max; %#ok<AGROW>
                ks0(i) = this.map(lbl).init; %#ok<AGROW>
            end
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end