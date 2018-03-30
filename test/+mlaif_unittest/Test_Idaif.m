classdef Test_Idaif < matlab.unittest.TestCase
	%% TEST_IDAIF 

	%  Usage:  >> results = run(mlaif_unittest.Test_Idaif)
 	%          >> result  = run(mlaif_unittest.Test_Idaif, 'test_dt')
 	%  See also:  file:///Applications/Developer/MATLAB_R2014b.app/help/matlab/matlab-unit-test-framework.html

	%  $Revision$
 	%  was created 17-Jan-2016 17:57:19
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlaif/test/+mlaif_unittest.
 	%% It was developed on Matlab 8.5.0.197613 (R2015a) for MACI64.
 	

	properties
 		registry
 		testObj
 	end

	methods (Test)
		function test_afun(this)
 			import mlaif.*;
 			this.assumeEqual(1,1);
 			this.verifyEqual(1,1);
 			this.assertEqual(1,1);
 		end
	end

 	methods (TestClassSetup)
		function setupIdaif(this)
 			import mlaif.*;
 			this.testObj_ = Idaif;
 		end
	end

 	methods (TestMethodSetup)
		function setupIdaifTest(this)
 			this.testObj = this.testObj_;
 		end
	end

	properties (Access = 'private')
 		testObj_
 	end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

