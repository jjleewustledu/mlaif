classdef AifPET < mlaif.PETDynamicsAbstract
	%% AIFPET imports, stores, interpolates, models, transforms arterial input functions for PET.
    %
	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.3.0.532 (R2014a) 
 	%  $Id$ 
    
    properties
        displayPETCurves = false
        kdefault = 0.01
        lambda = 0.9
        nativeMeasurements
        nativeTimes
        tracerConcentrations
    end
 	 
	properties (Dependent) 	
        aifMR
 	end 
    
    methods %% GET
        function a  = get.aifMR(this)
            assert(~isempty(this.aifMR_));
            a = this.aifMR_;
        end
    end
    
	methods (Static)
        function this = load(dcvName, perfLogName, varargin)
            %% LOAD
            %  Usage:   this = AifPET.load(crv_filename, perfusion4dfp_logfile[, parameter_name, parameter_value ...])
            
            this  = mlaif.AifPET(dcvName, perfLogName, varargin{:});
        end
    end 
    
    methods 
        function c              = CBF(this, k)
            if (~exist('k', 'var')); k = this.kdefault; end
            c = this.aifMR.CBF * k;
        end
        function c              = CBV(this, k)
            if (~exist('k', 'var')); k = this.kdefault; end
            c = this.CBF(k) * this.MTT;
        end
        function m              = MTT(this)
            m = 1 / this.aifMR.delta;
        end
        function [petaif,wbtac] = doublePass(this, params)
            %% DOUBLEPASS
            %  function to calculate the PET time courses from MRI and param data
            
            % PET AIF from MRI parameters
            petaif = params.ncnt * conv(this.aifMR.doublePass(params), this.mrToPetKernel(params));
            petaif = this.ensureSegment(petaif, 0, this.nTimes);

            % convolve PET residue curves to get WB pet signal            
            wbtac   = conv( ...
                          this.ensureLength(petaif, 2*this.nTimes-1), ...
                          this.aifToWholebrainKernel(params)) / params.fracaif;
            wbtac   = this.ensureSegment(wbtac, params.toffset, this.nTimes);
        end        
        function [petaif,wbtac] = doublePassFast(this, c1, c2, c3, c4, fracSS, fracaif, k, kappa, ncnt, t01, t02, toffset)
            %% DOUBLEPASS
            %  function to calculate the PET time courses from MRI and param data
            
            % PET AIF from MRI parameters
            petaif = ncnt * conv(this.aifMR.doublePassFast(fracSS, t01, t02), this.mrToPetKernelFast(kappa, c1, c2));
            petaif = this.ensureSegment(petaif, 0, 120);

            % convolve PET residue curves to get WB pet signal            
            wbtac   = conv( ...
                          this.ensureLength(petaif, 239), ...
                          this.aifToWholebrainKernelFast(k, c3, c4)) / fracaif;
            wbtac   = this.ensureSegment(wbtac, toffset, 120);
        end
        function [petaif,wbtac] = singlePass(this, params)
            %% SINGLEPASS
            %  function to calculate the PET time courses from MRI and param data
            
            % PET AIF from MRI parameters
            petaif_ = conv(this.aifMR.singlePass(params), this.mrToPetKernel(params));
            petaif  = this.ensureLength(petaif_, this.nTimes);

            % convolve PET residue curves to get WB pet signal            
            petaif_ = this.ensureLength(petaif_, this.nTimes+params.toffset);
            wbtac   = conv(petaif_, this.aifToWholebrainKernel(params));
            wbtac   = this.ensureSegment(wbtac, params.toffset, this.nTimes);

            petaif = params.ncnt*params.fracaif*petaif;
            wbtac  = params.ncnt*wbtac;
        end
        function d              = diracDelta(this)
            d    = zeros(1, this.nTimes);
            d(1) = 1;
        end
        function kernel         = mrToPetKernelPrevious(this, params)
            cbf    = this.CBF(params.k);
            cbv    = this.CBV(params.k);
            kernel = this.dt*cbf*(1/cbv - 1/this.lambda)*exp(-cbf*(1/cbv - 1/this.lambda)*this.timeInterpolants);
            kernel = kernel / max(kernel);
        end
        function kernel         = mrToPetKernel(this, params)
            arg    = -(1/params.kappa) * this.timeInterpolants;
            polyn  = 1 + params.c1 * arg + params.c2 * arg.^2;
            kernel = exp(arg) .* polyn;
            kernel = kernel / max(kernel);
        end
        function kernel         = mrToPetKernelFast(~, kappa, c1, c2)
            TIMES  = 0:119;
            arg    = -(1/kappa) * TIMES;
            polyn  = 1 + c1 * arg + c2 * arg.^2;
            kernel = exp(arg) .* polyn;
            kernel = kernel / max(kernel);
        end
        function kernel         = aifToWholebrainKernel(this, params)
            arg   = -this.CBF(params.k) * this.timeInterpolants;
            polyn = 1 + params.c3 * arg + params.c4 * arg.^2 + params.c5 * arg.^3 + params.c6 * arg.^4;
            kernel = exp(arg) .* polyn;
            kernel = kernel / max(kernel);
        end
        function kernel         = aifToWholebrainKernelFast(this, k, c3, c4)
            TIMES  = 0:119;
            arg    = -this.aifMR.CBF * k * TIMES;
            polyn = 1 + c3 * arg + c4 * arg.^2;
            kernel = exp(arg) .* polyn;
            kernel = kernel / max(kernel);
        end
    end
    
    %% PRIVATE
    
    properties (Access = 'private')
        aifMR_
        dcv_
    end
    
    methods (Static, Access = 'private')
        function f = ensureLength(f, len)
            %% ENSURELENGTH pads the end of a vector with the last vector value as needed to ensure a requested length
            %  padded_vec = this.ensureLength(vec, requested_length)
            
            lenf = length(f);
            len  = ceil(len);
            if (lenf < len)
                f(lenf+1:len) = f(end)*ones(1, len-lenf);
            end
            if (len < lenf)
                f = f(1:len);
            end
        end 
        function f = ensureSegment(f, start, len)
            start = ceil(start);
            len   = ceil(len);
            f     = mlaif.AifPET.ensureLength(f, len+start);
            try
                f     = f(start+1:start+len);
            catch ME
                error('mlaif:paramValuesOutOfBounds', 'AifPET.ensureSegment:  start->%i, len->%i', start, len);
            end
        end
    end
    
    methods (Access = 'private')
        function this = AifPET(dcvName, perfLogName, varargin)
 			%% AIFPET
 			%  Usage:   this = AifPET([parameter, parameter_values, ...])
            %                          ^ TimeInterpolants            
            
            this = this@mlaif.PETDynamicsAbstract(varargin{:});
            
            import mlaif.*;
            this.aifMR_ = AifMR.load(perfLogName);
            [dcvDir,dcvName] = fileparts(dcvName);
            pwd0 = pwd;
            cd(dcvDir);
            this.dcv_ = mlpet.CRV(dcvName);
            cd(pwd0);
            this.nativeMeasurements = this.clipFirstMeasurement(this.dcv_.counts');
            this.nativeTimes = this.dcv_.times'; 
            
            this.tracerConcentrations = this.clipFirstMeasurement(this.nativeMeasurements);
            if (~this.equivalent(this.timeInterpolants, this.nativeTimes))
                this.tracerConcentrations = this.pchip; end
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

