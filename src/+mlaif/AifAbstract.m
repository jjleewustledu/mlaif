classdef AifAbstract < mlaif.AifInterface 
	%% AIFABSTRACT manages interpolation of time
    %
	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.3.0.532 (R2014a) 
 	%  $Id$ 
 	 
	properties (Dependent)
        defaultTimeInterpolants           
        dt % of timeInterpolant; modal value if dt = f(t) or otherwise non-uniform        
        nTimes % of timeInterpolant
        timeInterpolants
    end

    methods %% GET/SET
        function d    = get.defaultTimeInterpolants(this)
            assert(~isempty(this.defaultTimeInterpolants_));
            d = this.defaultTimeInterpolants_;
        end
        function this = set.defaultTimeInterpolants(this, ti)
            assert(isnumeric(ti));
            assert(size(ti,2) > 1); % ensure nontrivial row vector
            this.defaultTimeInterpolants_ = ti;
        end
        function d    = get.dt(this)
            d = this.timeInterpolants(2:end) - this.timeInterpolants(1:end-1);
            d = mode(d);
        end
        function n    = get.nTimes(this)
            n = length(this.timeInterpolants);
        end
        function t    = get.timeInterpolants(this)
            assert(~isempty(this.timeInterpolants_));
            t = this.timeInterpolants_;
        end
    end
    
    %% PROTECTED
    
	properties (Access = 'protected') 		 
        defaultTimeInterpolants_ = 0:1:119
        timeInterpolants_
 	end 

	methods (Access = 'protected')
        function this = AifAbstract(varargin)            
 			%% AIFABSTRACT 
 			%  Usage:   this = AifAbstract([parameter, parameter_values, ...])
            %                               ^ TimeInterpolants

            p = inputParser;
            addParameter(p, 'TimeInterpolants', this.defaultTimeInterpolants, @isnumeric);
            parse(p, varargin{:});
            this.timeInterpolants_ = p.Results.TimeInterpolants;
        end
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

