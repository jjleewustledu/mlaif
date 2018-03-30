classdef AbstractVectorAifProblem < mlbayesian.AbstractVectorBayesianProblem
	%% ABSTRACTVECTORAIFPROBLEM   

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.3.0.532 (R2014a) 
 	%  $Id$ 
    
    properties (Abstract)
        baseTitle
        xLabel
        yLabel
    end
    
    properties (Constant)
        PLOT_ORI      = true
        PLOT_ESTIMATE = true
    end
    
    properties (Dependent)
        length
        timeInterpolants
    end
    
    methods %% GET
        function le = get.length(this)
            le = length(this.independentData);
        end
        function ti = get.timeInterpolants(this)
            ti = this.independentData;
        end
    end

    methods
        function this = AbstractVectorAifProblem(varargin)
 			%% ABSTRACTVECTORAIFPROBLEM 
 			%  Usage:  this = AbstractAifProblem([independent_data, dependent_data]) 
            
            this = this@mlbayesian.AbstractVectorBayesianProblem(varargin{:});
        end
        function this = runMcmc(this, paramsMap)
            %% RUNMCMC 
            %  Usage:  map = containers.Map
            %          map('param') = struct('fixed', 0, 'min', eps, 'mean', 1, 'max',  10)
            %          this = this.runMcmc(map)
            
            import mlbayesian.*;
            this.paramsManager = VectorBayesianParameters(paramsMap);            
            this.mcmc          = MCMC(this, this.dependentData, this.paramsManager);
            [~,~,this.mcmc]    = this.mcmc.runMcmc;            
                      
            this.plotOri; 
            this.plotEstimate            
        end   
        function        plotOri(this)
            if (~this.PLOT_ORI)
                return; end
            figure
            plot(this.timeInterpolants, this.dependentData, 'k', 'LineWidth', 2)
            title(this.baseTitle)
            xlabel(this.xLabel)
            ylabel(this.yLabel)
        end
        function        plotEstimate(this)
            if (~this.PLOT_ESTIMATE)
                return; end
            figure
            plot(this.timeInterpolants, this.dependentData, 'k', ...
                 this.timeInterpolants, this.estimateData,  'k:', 'LineWidth', 2);
            title(sprintf('%s and Bayesian estimate', this.baseTitle));
            xlabel(this.xLabel)
            ylabel(this.yLabel)
        end
    end
    
    %% PROTECTED
    
    methods (Access = 'protected')
        function [c1,c2] = twinBolusPassagesFast(this, a, b, t0, tauRe)
            idx_t0 = round(t0) + 1;
            idx_t1 = round(t0 + tauRe) + 1;
            times  = this.timeInterpolants;
            normg  = gamma(a + 1) / b^(a + 1);
            c0     = times.^a .* exp(-b * times) / normg; 
            c1     = zeros(1, length(times));
            c2     = zeros(1, length(times));              
            c1(idx_t0:end) = c0(1:end-idx_t0+1);
            c2(idx_t1:end) = c0(1:end-idx_t1+1);
        end           
        function  c      = bolusSteadyStateFast( this, g, t0)
            idx_t0 = floor(t0) + 1;
            times  = this.timeInterpolants;
            c0     = (1 - exp(-g * times));   
            c      = zeros(1, length(times));         
            c(idx_t0:end) = c0(1:end-idx_t0+1);
        end
        
        function v = betaVariateFast(     this, x, y, t0)            
            %% GAMMAVARIATEFAST
            %  v(t,x,y) = \frac{ (t - t0)^{x - 1} }{ Beta(x, y) (1 + t - t0)^{x + y} },
            %             Re{x} > 0, Re{y} > 0
            
            idx_t0 = floor(t0) + 1;
            times  = this.timeInterpolants;
            v0     = times.^(x - 1) ./ (1 + times).^(x + y) / beta(x, y);
            v      = zeros(1, length(times));    
            v(idx_t0:end) = v0(1:end-idx_t0+1);
        end
        function v = gammaVariatePolyFast(this, a, b, c1, c2, c3, c4, t0)
            %% GAMMAVARIATEPOLYFAST
            %  arg(b,t0)   =      b (t - t0)
            %  v(t,a,b,t0) = \frac{ (t - t0)^a exp(-arg) [1 + c1 arg + c2 arg^2 + c3 arg^3 + c4 arg4]}{ Gamma(a + 1) b^{a + 1} }, 
            %                a > -1 to avoid poles, t0 >= 0
            
            idx_t0 = round(t0) + 1;
            times  = this.timeInterpolants;
            normg  = gamma(a + 1)/b^(a + 1);
            arg    = b * times;
            poly   = 1 + c1 * arg + c2 * arg.^2 + c3 * arg.^3 + c4 * arg.^4;
            v0     = times.^a .* exp(-arg) .* poly / normg; 
            v      = zeros(1, length(times)); 
            v(idx_t0:end) = v0(1:end-idx_t0+1);
        end
        function v = gammaVariateFast(    this, a, b, t0)
            %% GAMMAVARIATEFAST
            %  v(t,a,b,t0) = \frac{ (t - t0)^a exp(-b (t - t0)) }{ Gamma(a + 1) b^{a + 1} }, 
            %                a > -1 to avoid poles, t0 >= 0
            
            idx_t0 = floor(t0) + 1;
            times  = this.timeInterpolants;
            normg  = gamma(a + 1) / b^(a + 1);
            v0     = times.^a .* exp(-b * times)/normg;      
            v      = zeros(1, length(times));    
            v(idx_t0:end) = v0(1:end-idx_t0+1);
        end
        function c = expPolyFast(         this, c1, c2, c3, c4, d, t0)
            %% EXPPOLYFAST
            %  c(t,d,t0) = exp(-arg) [1 + c1 arg + c2 arg^2 + c3 arg^3  + c4 arg^4], arg = d (t - t0)
            
            idx_t0 = floor(t0) + 1;
            times  = this.timeInterpolants;
            arg    = d * times;
            c0     = exp(-arg) .* (1 + c1 * arg + c2 * arg.^2 + c3 * arg.^3 + c4 * arg.^4);
            c      = zeros(1, length(times));
            c(idx_t0:end) = c0(1:end+1-idx_t0);
        end   
        function c = expFast(             this, d,  t0)
            %% EXPFAST
            %  c(t,d,t0) = exp(-d (t - t0))
            
            idx_t0 = floor(t0) + 1;
            times  = this.timeInterpolants;
            c0     = exp(-d * times);
            c      = zeros(1, length(times));
            c(idx_t0:end) = c0(1:end-idx_t0+1);
        end    
    end
    
	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

