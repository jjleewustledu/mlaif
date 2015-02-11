classdef Test_BayesianPET < matlab.unittest.TestCase 
	%% TEST_BAYESIANPET  

	%  Usage:  >> results = run(mlaif_unittest.Test_BayesianPET)
 	%          >> result  = run(mlaif_unittest.Test_BayesianPET, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014a.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.4.0.150421 (R2014b) 
 	%  $Id$ 
 	 

	properties 
        bayesianPET
 	end 

	methods (Test)
 		function test_estimateParameters(this)  			
            this.assertNotEmpty(this.bayesianPET.estimateParameters);
 		end 
    end 
    
    methods (TestClassSetup)
        function setupBayesianPET(this)
            this.bayesianPET = mlaif.BayesianPET(fullfile('/Volumes/PassportStudio2/cvl/np755/mm01-007_p7267_2008jun16'), 'p7267ho1');
        end
    end
    
    methods (TestClassTeardown)
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

