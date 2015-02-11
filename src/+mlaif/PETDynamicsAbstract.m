classdef PETDynamicsAbstract < mlaif.AifAbstract
	%% PETDYNAMICSABSTRACT defines tracerConcentrations, spline so as to be useful for most
    %  PET time-series.
	%
	%  $Revision$ 
	%  was created $Date$ 
	%  by $Author$, 
	%  last modified $LastChangedDate$ 
	%  and checked into repository $URL$, 
	%  developed on Matlab 8.3.0.532 (R2014a) 
	%  $Id$ 
	%  N.B. classdef (Sealed, Hidden, InferiorClasses = {?class1,?class2}, ConstructOnLoad) 

    methods   
        function m    = clipFirstMeasurement(~, m)
            baseline = mean(m(3:10));
            if (abs(m(1)/baseline) > 10)
                m(1) = m(3);
                m(2) = m(3);
            end
        end     
 		function yy = spline(this)
            if (this.timeInterpolants(end) > this.nativeTimes(end))
                error('mlaif:notImplemented', ...
                      'requested time-point of cubic spline %g exceeds longest time-point of native data %g', ...
                      this.timeInterpolants(end), this.nativeTimes(end));
            end
            Y          = zeros(1, length(this.nativeMeasurements)+2);
            Y(2:end-1) = this.nativeMeasurements(1:end);
            yy         = spline(this.nativeTimes, Y, this.timeInterpolants);
        end
    end
    
    %% PROTECTED
    
    methods (Static, Access = 'protected')        
        function tf = equivalent(a, b)
            tf = true;
            assert(length(a) == numel(a));
            assert(length(b) == numel(b));
            if (any(size(a) ~= size(b)))
                tf = false; return; end
            if (any(isnan(a)) || any(isnan(b)))
                tf = false; return; end
            if (any(a ~= b))                
                tf = false; return; end
        end
    end
    
    methods (Access = 'protected')   
 		function this = PETDynamicsAbstract(varargin)
			%% PETDYNAMICSABSTRACT 
 			%  Usage:   this = PETDynamicsAbstract([parameter, parameter_values, ...])
            %                                      ^ TimeInterpolants
            
            this = this@mlaif.AifAbstract(varargin{:});
		end      
    end
    
	%  Created with NewClassStrategy by John J. Lee, after newfcn by Frank Gonzalez-Morphy 
end

