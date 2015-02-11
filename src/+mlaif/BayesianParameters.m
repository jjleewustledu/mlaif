classdef BayesianParameters  
	%% BAYESIANPARAMETERS   

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.3.0.532 (R2014a) 
 	%  $Id$  	 

    properties
        max
        min
        fixed
        fixedValue
        mean
        std
        nTimes
 	end 

	properties (Dependent)
        length
    end
    
    methods %% GET/SET
        function len = get.length(this)
            len = length(this.mean);
        end
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

