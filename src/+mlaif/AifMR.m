classdef AifMR < mlaif.AifAbstract 
	%% AIFMR imports, stores, interpolates, models, transforms arterial input functions for MR.
    %  Data from this class are entirely obtained from Bayesian parameters estimated from native bolus-tracking data.
    %
	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.3.0.532 (R2014a) 
 	%  $Id$ 
 	 
    properties
        k = 6 % empirical estimate  
        nativeMeasurements
        nativeTimes 
        tracerConcentrations
    end
    
	properties (Dependent)        
        S0
        CBF
        CBV
        t0
        alpha
        beta
        delta
        gamma
        eps
    end 
    
    methods (Static)
        function this   = load(fname, varargin)
            %% LOAD
            %  Usage:   this = AifMR.load(perfusion4dfp_logfile[, parameter_name, parameter_value])
            
            this = mlaif.AifMR(fname, varargin{:});
        end 
    end
    
    methods %% GET
        function x = get.S0(this)
            x = this.perfusionIO_.S0;
        end
        function x = get.CBF(this)
            x = this.perfusionIO_.CBF;
        end
        function x = get.CBV(this)
            x = this.CBF/this.delta;
        end
        function x = get.t0(this)
            x = this.perfusionIO_.t0;
        end
        function x = get.alpha(this)
            x = this.perfusionIO_.alpha;
        end
        function x = get.beta(this)
            x = this.perfusionIO_.beta;
        end
        function x = get.delta(this)
            x = this.perfusionIO_.delta;
        end
        function x = get.gamma(this)
            x = this.perfusionIO_.gamma;
        end
        function x = get.eps(this)
            x = this.perfusionIO_.eps;
        end
        
    end
    
    methods
        function aif    = doublePass(this, params)
            %% DOUBLEPASS describe a double-pass including recirculation
            
            fss = params.fracSS;
            
            % MRI AIF main component
            ca1 = (1 - fss)*this.eps*this.bolusPassage(this.alpha, this.beta, params.t01);
            
            % recirculation part
            ca2 = (1 - fss)*(1 - this.eps)*this.bolusPassage(this.alpha, this.beta, params.t01 + params.t02);
            
            % steady state piece
            ca3 = max(ca1)*fss*this.bolusSteadyState(this.gamma, params.t01);

            % put all the pieces together
            aif = ca1 + ca2 + ca3;
        end  
        function aif    = doublePassFast(this, fracSS, t01, t02)
            %% DOUBLEPASSFAST describe a double-pass including recirculation
            
            [bp1,bp2] = this.twinBolusPassagesFast(this.alpha, this.beta, t01, t01+t02);
            
            % MRI AIF main component
            ca1 = (1 - fracSS)*this.eps*bp1;
            
            % recirculation part
            ca2 = (1 - fracSS)*(1 - this.eps)*bp2;
            
            % steady state piece
            ca3 = max(ca1)*fracSS*this.bolusSteadyStateFast(this.gamma, t01);

            % put all the pieces together
            aif = ca1 + ca2 + ca3;
        end  
        function aif    = singlePass(this, params)
            %% SINGLEPASS simplifies calc_mriaif to describe a single-pass without recirculation
            
            fss = params.fracSS;
            
            % MRI AIF main component
            ca1 = (1 - fss)*this.bolusPassage(this.alpha, this.beta, params.t01);

            % steady state piece
            ca3 = max(ca1)*fss*this.bolusSteadyState(this.gamma, params.t01);

            % put all the pieces together
            aif = ca1 + ca3;
        end    
        function aif    = singlePassFast(this, fracSS, t01)
            %% DOUBLEPASSFAST describe a double-pass including recirculation
            
            % MRI AIF main component
            ca1 = (1 - fracSS)*this.bolusPassagesFast(this.alpha, this.beta, t01);
            
            % steady state piece
            ca3 = max(ca1)*fracSS*this.bolusSteadyStateFast(this.gamma, t01);

            % put all the pieces together
            aif = ca1 + ca3;
        end  
        function conc   = bolusPassage(this, alph, bet, t0) 
            times = this.timeInterpolants;          
            conc  = zeros(1, length(times));
            norm  = gamma(alph + 1)/bet^(alph + 1);
            curv  = this.dt*(times.^alph).*exp(-bet*times)/norm;
            index0 = 0;
            for i = 1:length(times)
                if times(i) < t0
                    conc(i) = 0;
                    index0  = i;
                else
                    conc(i) = curv(i - index0);
                end
            end
        end
        function conc   = gammaVariateFast(~, alph, bet, t0)
            idx_t0 = round(t0) + 1;
            times  = 0:1:119;
            conc   = zeros(1, 120);
            norm   = gamma(alph + 1)/bet^(alph + 1);
            curv   = (times.^alph).*exp(-bet*times)/norm;               
            conc(idx_t0:end) = curv(1:end+1-idx_t0);
        end
        function [conc,conc1] = twinBolusPassagesFast(~, alph, bet, t0, t1)
            idx_t0 = round(t0) + 1;
            idx_t1 = round(t1) + 1;
            times  = 0:1:119;
            conc   = zeros(1, 120);
            conc1  = zeros(1, 120);
            norm   = gamma(alph + 1)/bet^(alph + 1);
            curv   = (times.^alph).*exp(-bet*times)/norm;               
            conc( idx_t0:end) = curv(1:end+1-idx_t0);
            conc1(idx_t1:end) = curv(1:end+1-idx_t1);
        end
        function  conc        = bolusPassagesFast(~, alph, bet, t0)
            idx_t0 = round(t0) + 1;
            times  = 0:1:119;
            conc   = zeros(1, 120);
            norm   = gamma(alph + 1)/bet^(alph + 1);
            curv   = (times.^alph).*exp(-bet*times)/norm;               
            conc(idx_t0:end) = curv(1:end+1-idx_t0);
        end
        function conc   = bolusSteadyState(this, gam, t0)
            times = this.timeInterpolants;
            conc  = zeros(1, length(times));
            curv  = this.dt*(1 - exp(-gam*times));
            for i = 1:length(times)
                if times(i) < t0
                    conc(i) = 0;
                    index0  = i;
                else
                    conc(i) = curv(i - index0);
                end
            end
        end
        function conc   = bolusSteadyStateFast(~, gam, t0)
            idx_t0 = round(t0) + 1;
            times  = 0:1:119;
            conc   = zeros(1, 120);
            curv   = (1 - exp(-gam*times));            
            conc(idx_t0:end) = curv(1:end+1-idx_t0);
        end
        function aif    = singlePassDebug(this)
            %% SINGLEPASS simplifies calc_mriaif to describe a single-pass without recirculation            
            
            ti    = this.timeInterpolants;
            t1    = 12.418047985835880;
            t2    = 28.602017457330973;
            a1    = 7.833932447361581;
            a2    = 6.551300832188158;
            b1    = 0.594878893044038;
            b2    = 0.218221163899827;
            q2    = 0.306217434449070;
            q3    = 0.0015;
            
            % MRI AIF main component
            ca1   = zeros(1, length(ti));
            n1    = gamma(a1+1.0)/b1^(a1+1.0);
            temp1 = this.dt*(ti.^a1).*exp(-b1*ti)/n1;
            for i = 1:length(ti)
                if ti(i)<t1
                    ca1(i) = 0.0;
                    i1 = i;
                else
                    ca1(i) = temp1(i-i1);
                end
            end
            
            % recirculation part
            ca2   = zeros(1, length(ti));
            n2    = gamma(a2+1.0)/b2^(a2+1.0);
            temp2 = q2*this.dt*(ti.^a2).*exp(-b2*ti)/n2;
            for i = 1:length(ti)
                if ti(i)<t2
                    ca2(i) = 0.0;
                    i2 = i;
                else
                    ca2(i) = temp2(i-i2);
                end
            end

            % steady state piece
            sumca = zeros(1, length(ti));
            sumca(1) = ca1(1);
            for i = 2:length(ti)
                sumca(i) = sumca(i-1) + q3*ca1(i); % compare to q3 below
            end

            % put all the pieces together
            aif = ca1 + ca2 + sumca;
        end
    end

    %% PRIVATE
    
    properties (Access = 'private')       
        perfusionIO_ % logfile from Bayesian perfusion4dfp
    end
    
	methods (Access = 'private')
 		function this   = AifMR(fname, varargin) 
 			%% AIFMR 
 			%  Usage:   this = AifMR([parameter, parameter_values, ...])
            %                         ^ TimeInterpolants

            this = this@mlaif.AifAbstract(varargin{:});   
            
            this.perfusionIO_ = mlperfusion.PerfusionIO.load(fname);
            
            this.tracerConcentrations = this.doublePass( ...
                                        struct('t01', this.t0, 't02', 10, 'fracSS', 0.17));
            this.nativeMeasurements = zeros(1, this.nTimes);
            for n = 1:this.nTimes
                this.nativeMeasurements(n) = this.S0 * exp(-this.k*this.dt*this.tracerConcentrations(n));
            end                
            this.nativeTimes = this.timeInterpolants;
                        
 		end         
 	end 

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

 