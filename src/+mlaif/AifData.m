classdef AifData < handle
    %% AIFDATA is a data class gathering information needed by preprocessors and models.
    %  To do:  convert from singleton to simple handle.
    %  
    %  Created 19-Apr-2022 14:30:05 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlaif/src/+mlaif.
    %  Developed on Matlab 9.12.0.1884302 (R2022a) for MACI64.  Copyright 2022 John J. Lee.
    
    methods (Static)
        function this = instance(varargin)
            %% INSTANCE
            %  @param optional qualifier is char \in {'initialize' ''}
            
            ip = inputParser;
            addOptional(ip, 'qualifier', '', @ischar)
            parse(ip, varargin{:})
            
            persistent uniqueInstance
            if (strcmp(ip.Results.qualifier, 'initialize'))
                uniqueInstance = [];
            end          
            if (isempty(uniqueInstance))
                this = mlaif.AifData();
                uniqueInstance = this;
            else
                this = uniqueInstance;
            end
        end
    end 

    properties (Dependent)
        Ddatetime0 % seconds.  scannerDev.datetime0 - arterialDev.datetime0, typically < 0 for automated arterial sampling.
        normalizationFactor % recovery coefficient for IDIF; empirical rescaling for debugging measurements of AIF.
        stableToInterpolation % logical.  Indicates that makima, pchip and related interpolators are safe for use.
        T % seconds.  The duration at the start of artery_interpolated used by models but not recorded by scanner time-frames.
        tArterialForced % prefer TwiliteDevice.t0_forced, RadialArteryLee2021Model.t0_forced.
        tBuffer % seconds.  artery_interpolated(tBuffer + 1) corresponds to the first scanner time-frame.
    end

    methods

        %% GET/SET

        function g = get.Ddatetime0(this)
            g = this.Ddatetime0_;
        end
        function g = get.normalizationFactor(this)
            g = this.normalizationFactor_;
        end
        function g = get.stableToInterpolation(this)
            g = this.stableToInterpolation_;
        end
        function g = get.T(this)
            g = this.T_;
        end
        function g = get.tArterialForced(this)
            g = this.tArterialForced_;
        end
        function g = get.tBuffer(this)
            g = max(0, -this.Ddatetime0) + this.T;
        end

        function set.Ddatetime0(this, s)
            assert(isscalar(s) || isduration(s))
            if isduration(s)
                s = seconds(s);
            end
            this.Ddatetime0_ = s;
        end
        function set.stableToInterpolation(this, s)
            assert(islogical(s))
            this.stableToInterpolation_ = s;
        end
        function set.normalizationFactor(this, s)
            assert(isscalar(s))
            this.normalizationFactor_ = s;
        end
        function set.T(this, s)
            assert(isscalar(s))
            this.T_ = s;
        end
        function set.tArterialForced(this, s)
            assert(isscalar(s));
            this.tArterialForced_ = s;
        end
    end
        
    %% PRIVATE

    properties (Access = private)
        Ddatetime0_
        normalizationFactor_
        stableToInterpolation_
        T_
        tArterialForced_
    end

    methods (Access = private)

        function this = AifData(varargin)
            ip = inputParser;
            addParameter(ip, 'Ddatetime0', 0, @isscalar);
            addParameter(ip, 'normalizationFactor', 1, @isscalar);
            addParameter(ip, 'stableToInterpolation', false, @islogical)
            addParameter(ip, 'T', 0, @isscalar);
            addParameter(ip, 'tArterialForced', 0, @isscalar);
            parse(ip, varargin{:});
            ipr = ip.Results;

            this.Ddatetime0_ = ipr.Ddatetime0;
            this.normalizationFactor_ = ipr.normalizationFactor;
            this.stableToInterpolation_ = ipr.stableToInterpolation;
            this.T_ = ipr.T;
            this.tArterialForced_ = ipr.tArterialForced;
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
