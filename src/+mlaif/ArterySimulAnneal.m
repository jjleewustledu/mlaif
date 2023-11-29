classdef ArterySimulAnneal < mlio.AbstractIO
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

        context
        Measurement          % external data as table
        map                  % 
        model                % mlaif.ArteryModel
        sigma0 = 0.05        % fraction of Measurement < 1
    end
    
	properties (Dependent) 
        kernel  
        ks
        ks_names
        model_kind
        product
        tracer
    end

	methods %% GET
        function g = get.kernel(this)
            g = this.model.kernel;
        end
        function g = get.ks(this)
            g = this.product_.ks;
        end
        function g = get.ks_names(this)
            g = this.model.ks_names;
        end
        function g = get.model_kind(this)
            g = this.model.model_kind;
        end
        function g = get.product(this)
            g = this.product_;
        end
        function g = get.tracer(this)
            g = this.model.tracer;
        end
    end

    methods
        function fprintfModel(this)
            fprintf('ArterySimulAnneal:\n');  
            fprintf('%s %s:\n', this.tracer, this.model_kind); 
            for ky = 1:length(this.ks)
                fprintf('\t%s = %g\n', this.ks_names{ky}, this.ks(ky));
            end 
            fprintf('\tloss = %g\n', this.loss())
            keys = natsort(this.map.keys);
            for ky = 1:length(this.ks)
                fprintf('\tmap(''%s'') => %s\n', this.ks_names{ky}, ...
                    join(struct2str(this.map(keys{ky}), orientation='horz')));
            end
        end
        function Q = loss(this)
            Q = this.product_.loss;
        end
        function h = plot(this, opts)
            arguments
                this mlaif.ArterySimulAnneal
                opts.showKernel logical = true
                opts.xlim double = [-10 200]
                opts.ylim double = []
                opts.zoom double = []
            end
            this.zoom = opts.zoom;
            timesMid = asrow(this.Measurement.timesMid);
            measurement = asrow(this.Measurement.activityDensity);
            N = length(measurement);
            M0 = max(measurement);                       
            samp = this.model.sampled(M0, this.ks, N, this.kernel, this.tracer, this.model_kind);
            deconvolved = this.model.deconvolved(M0, this.ks, N, this.kernel, this.tracer, this.model_kind);
            
            if isempty(this.zoom)
                % show kernel as 1/2 amplitude of deconvolved
                this.zoom = max(deconvolved)/max(this.kernel)/2; 
            end
            if this.zoom ~= 1
                leg_kern = sprintf('kernel x%g', this.zoom);
            else
                leg_kern = 'kernel';
            end
 
            h = figure;
            if opts.showKernel && ~isscalar(this.kernel)
                hold('on')
                plot(timesMid, measurement, 'o', 'MarkerEdgeColor', "#0072BD")
                plot(timesMid, samp, '--', 'Color', "#A2142F", 'LineWidth', 2)
                plot(timesMid, deconvolved, '-', 'Color', "#0072BD", 'LineWidth', 2)
                plot(timesMid, this.zoom*this.kernel(1:N), '--', 'Color', "#EDB120", 'LineWidth', 2)
                legend({'measured', 'estimated', 'deconvolved', leg_kern}, 'FontSize', 12)
                hold('off')
            else
                hold('on')
                plot(timesMid, measurement, 'o', 'MarkerEdgeColor', "#0072BD")
                plot(timesMid, samp, '--', 'Color', "#A2142F", 'LineWidth', 2)
                plot(timesMid, deconvolved, '-', 'Color', "#0072BD", 'LineWidth', 2)
                legend({'measured', 'estimated', 'deconvolved'}, 'FontSize', 12)
                hold('off')
            end
            if ~isempty(opts.xlim); xlim(opts.xlim); end
            if ~isempty(opts.ylim); ylim(opts.ylim); end
            xlabel('times / s', FontSize=14, FontWeight='bold')
            ylabel('activity / (Bq/mL)', FontSize=14, FontWeight='bold')
            annotation('textbox', [.25 .5 .3 .3], 'String', sprintfModel(this), 'FitBoxToText', 'on', 'FontSize', 10, 'LineStyle', 'none')
            title(clientname(false, 2), FontSize=14)
            set(h, position=[100,100,1000,618])
        end
        function save(this)
            save([this.fileprefix '.mat'], this);
        end
        function saveas(this, fn)
            save(fn, this);
        end
        function this = solve(this, loss_function)
            %% Args:
            %      this mlaif.ArterySimulAnneal
            %      loss_function function_handle
            
            arguments
                this mlaif.ArterySimulAnneal
                loss_function function_handle
            end
            
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
 			[ks_,loss,exitflag,output] = simulannealbnd( ...
                @(ks__) loss_function( ...
                       ks__, this.kernel, this.tracer, this.model_kind, this.Measurement), ...
                this.ks0, this.ks_lower, this.ks_upper, options); 
            
            this.product_ = struct('ks0', this.ks0, 'ks', ks_, 'loss', loss, 'exitflag', exitflag, 'output', output); 
            if ~this.quiet
                fprintfModel(this)
            end
            if this.visualize
                plot(this)
            end
        end 
        function s = sprintfModel(this)
            s = sprintf('ArterySimulAnnel:\n');
            s = [s sprintf('%s %s:\n', this.tracer, this.model_kind)];
            for ky = 1:length(this.ks)
                s = [s sprintf('\t%s = %g\n', this.ks_names{ky}, this.ks(ky))]; %#ok<AGROW>
            end
            s = [s sprintf('\tloss = %g\n', this.loss())];
            keys = natsort(this.map.keys);
            for ky = 1:length(this.ks)
                s = [s sprintf('\tmap(''%s'') => %s\n', this.ks_names{ky}, ...
                    join(struct2str(this.map(keys{ky}), orientation='horz')))]; %#ok<AGROW>
            end
        end

        function this = ArterySimulAnneal(opts)
            %%  Args:
            %   context
            %   opts.sigma0 double {mustBeScalarOrEmpty} = []
            
            arguments
                opts.context
                opts.sigma0 double {mustBeScalarOrEmpty} = []
            end
            if isempty(opts.sigma0)
                opts.sigma0 = this.sigma0;
            end
            
            this.context = opts.context;
            this.model = this.context.model;               % copy objects for speed            
            this.map = this.model.map;                     %
            this.Measurement = this.context.Measurement;   %
            this.sigma0 = opts.sigma0;

            [this.ks_lower,this.ks_upper,this.ks0] = remapper(this);
 		end
    end
    
    %% PROTECTED
    
    properties (Access = protected)
        product_
    end
    
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
