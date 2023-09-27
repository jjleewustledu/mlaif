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
                dt = this.json_.start_times(2) - this.json_.start_times(1);
                timesF = cumsum(this.json_.taus);
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
        function q = deconv(this)
        end
        function q = deconvBayes(this)
        end
        function k = kernel(this)
            k = zeros(size(this.timeInterpolants));
            Nk = this.json_.taus(1);
            k(1:Nk) = 1/Nk;
        end
        function plot_deconv(this)
        end

    end

    methods (Static)
        function this = create(Measurement, opts)
            arguments
                Measurement {mustBeNonempty}
                opts.json struct = struct([])
                opts.tracer {mustBeText}
            end

            this = mlaif.Idif();
            this.Measurement_ = mlfourd.ImagingContext2(Measurement);
            if ~isempty(this.Measurement_.json_metadata)
                this.json_ = this.Measurement_.json_metadata;
            end
            if ~isempty(opts.json)
                this.json_ = opts.json;
            end
            if contains(this.Measurement_.fileprefix, "_trc-")
                re = regexp(this.Measurement_.fileprefix, "\S+_trc-(?<trc>[a-zA-Z0-9]+)_\S+", "names");
                this.tracer_ = re.trc;
            end
            if ~isempty(opts.tracer)
                this.tracer_ = opts.tracer;
            end
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
