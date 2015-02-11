classdef Test_AifMR < matlab.unittest.TestCase 
	%% TEST_AIFMR ...
    %  Usage:  >> results = run(mlaif_unittest.Test_AifMR)
	%          >> result  = run(mlaif_unittest.Test_AifMR, 'test_dt')
	%  See also:  file:///Applications/Developer/MATLAB_R2014a.app/help/matlab/matlab-unit-test-framework.html
    %
	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.3.0.532 (R2014a) 
 	%  $Id$ 

    properties
        testData
    end
    
	methods (Test)
        function test_dt(this)
            this.verifyEqual(this.testData.aifMR.dt, 1);
        end 
        function test_nTimes(this)
            this.verifyEqual(this.testData.aifMR.nTimes, 120);
        end 
        function test_timeInterpolants(this)
            this.verifyEqual(this.testData.aifMR.timeInterpolants, 0:119);
        end 
        
        function test_nativeMeasurements(this)
            nm = this.testData.aifMR.nativeMeasurements;
            this.verifyNumElements(nm, 120);
            this.verifyEqual(nm(1),   589.488);
            this.verifyEqual(nm(2),   589.488);
            this.verifyEqual(nm(3),   589.488);
            this.verifyEqual(nm(10),  589.488);
            this.verifyEqual(nm(120), 464.8423404537158, 'AbsTol', 1e-6);
        end 
        function test_nativeTimes(this)
            this.verifyEqual(this.testData.aifMR.nativeTimes, 0:119);    
        end 
        function test_tracerConcentrations(this)
            tc = this.testData.aifMR.tracerConcentrations;
            this.verifyNumElements(tc, 120);
            this.verifyEqual(tc(1),   0);
            this.verifyEqual(tc(2),   0);
            this.verifyEqual(tc(3),   0);
            this.verifyEqual(tc(10),  0);
            this.verifyEqual(tc(120), 0.039592678037980, 'AbsTol', 1e-6);
        end 
        
        function test_S0(this)
            this.verifyEqual(this.testData.aifMR.S0, 589.488);
        end 
        function test_CBF(this)
            this.verifyEqual(this.testData.aifMR.CBF, 1.885);
        end
        function test_t0(this)
            this.verifyEqual(this.testData.aifMR.t0, 8.152);
        end
        function test_alpha(this)
            this.verifyEqual(this.testData.aifMR.alpha, 6.031);
        end
        function test_beta(this)
            this.verifyEqual(this.testData.aifMR.beta, 1.866);
        end
        function test_delta(this)
            this.verifyEqual(this.testData.aifMR.delta, 0.725);
        end
        function test_gamma(this)
            this.verifyEqual(this.testData.aifMR.gamma, 1.702);
        end
        function test_eps(this)
            this.verifyEqual(this.testData.aifMR.eps, 0.954);
        end
        
        function test_gammaVariate(this)
            aifMR = this.testData.aifMR;
            alph = 1:2.718:10;
            bet  = 1:3.142:10;
            t0   = 0:10:60;
            for a = 1:length(alph)
                for b = 1:length(bet)
                    for t = 1:length(t0)
                        this.verifyEqual(aifMR.bolusPassage(    alph(a), bet(b), t0(t)), ...
                                         aifMR.gammaVariateFast(alph(a), bet(b), t0(t)), 'AbsTol', 1e-6);
                    end
                end
            end
        end
        function test_bolusSteadyStateFast(this)
            aifMR = this.testData.aifMR;
            alph = linspace(5,   9,   6);
            bet  = linspace(0.2, 0.8, 6);
            t0   = linspace(0,  60,   6);
            for a = 1:length(alph)
                for b = 1:length(bet)
                    for t = 1:length(t0)
                        this.verifyEqual(aifMR.bolusPassage(    alph(a), bet(b), t0(t)), ...
                                         aifMR.gammaVariateFast(alph(a), bet(b), t0(t)), 'AbsTol', 1e-6);
                    end
                end
            end
        end
        function test_singlePass(this)
            params.ncnt = 1e7;
            params.fracaif = 1;
            params.fracSS = 0;
            params.toffset = 10;
            params.t01 = 10;
            params.kappa = 1;
            params.c1 = 1;
            params.c2 = 1;
            params.c3 = 1;
            params.c4 = 1;  
            params.c5 = 1;
            params.c6 = 1;          
            params.k = 1;
            spp = this.testData.aifMR.singlePass(params);
            this.verifyNumElements(spp, 120);
            this.verifyEqual(spp(1),   0);
            this.verifyEqual(spp(2),   0);
            this.verifyEqual(spp(3),   0);
            this.verifyEqual(spp(10),  0);
            this.verifyEqual(spp(15), 0.2580189, 'AbsTol', 1e-6);
        end
    end
    
    methods (TestClassSetup)
        function addPathAndAifMR(this)
            this.testData.origPath = pwd;
            this.testData.workPath = fullfile( ...
                getenv('HOME'), ...
                'Local/src/mlcvl/mlaif/test/+mlaif_unittest');
            cd(this.testData.workPath);
            this.testData.aifMR = mlaif.AifMR.load( ...
                fullfile(this.testData.workPath, 'perfusion_4dfp.log'));
        end
    end
    
    methods (TestClassTeardown)
        function removePath(this)
            cd(this.testData.origPath);
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

