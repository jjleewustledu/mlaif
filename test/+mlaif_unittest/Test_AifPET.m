classdef Test_AifPET < matlab.unittest.TestCase  
	%% TEST_AIFPET  
    %  Usage:  >> results = run(mlaif_unittest.Test_AifPET)
	%          >> result  = run(mlaif_unittest.Test_AifPET, 'test_dt')
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
            this.verifyEqual(this.testData.aifPET.dt, 1);
        end 
        function test_nTimes(this)
            this.verifyEqual(this.testData.aifPET.nTimes, 120);
        end 
        function test_timeInterpolants(this)
            this.verifyEqual(this.testData.aifPET.timeInterpolants, 0:119);
        end 
        
        function test_nativeMeasurements(this)
            nm = this.testData.aifPET.nativeMeasurements;
            this.verifyNumElements(nm, 120);
            this.verifyEqual(nm(1),   14);
            this.verifyEqual(nm(2),   14);
            this.verifyEqual(nm(3),   48);
            this.verifyEqual(nm(10),  83);
            this.verifyEqual(nm(120), 2013);
        end 
        function test_nativeTimes(this)
            this.verifyEqual(this.testData.aifPET.nativeTimes, 1:120);    
        end 
        function test_tracerConcentrations(this)            
            tc = this.testData.aifPET.tracerConcentrations;
            this.verifyNumElements(tc, 120);
            this.verifyEqual(tc(1),  -29.497665284637662);
            this.verifyEqual(tc(2),   14);
            this.verifyEqual(tc(3),   14);
            this.verifyEqual(tc(10),  83);
            this.verifyEqual(tc(120), 1802);
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
            spp = this.testData.aifPET.singlePass(params);
            this.verifyNumElements(spp, 120);            
            this.verifyEqual(spp(1),   0);
            this.verifyEqual(spp(2),   0);
            this.verifyEqual(spp(3),   0);
            this.verifyEqual(spp(10),  0);
            this.verifyEqual(spp(20), 1451842, 'AbsTol', 1);
        end
    end
	
    methods (TestClassSetup)
        function addPathAndAifPET(this)
            this.testData.origPath = pwd;
            this.testData.workPath = fullfile( ...
                getenv('HOME'), ...
                'MATLAB-Drive/mlaif/test/+mlaif_unittest');
            cd(this.testData.workPath);
            this.testData.aifPET = mlaif.AifPET.load( ...
                fullfile(this.testData.workPath, 'p7267ho1.crv'), ...
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

