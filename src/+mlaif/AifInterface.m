classdef AifInterface  
	%% AIFINTERFACE specifies how to import, store, interpolate, model, transform arterial input functions 
    %
	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.3.0.532 (R2014a) 
 	%  $Id$  	 

	properties (Abstract)
        defaultTimeInterpolants
        dt % interpolant
        nTimes % interpolant        
        timeInterpolants
        
        nativeMeasurements
        nativeTimes
        tracerConcentrations
    end 

    methods (Abstract, Static)
        load % factory method; subclass ctors should be private
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

