classdef Test_MipIdif < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 21-Nov-2023 20:52:08 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlaif/test/+mlaif_unittest.
    %  Developed on Matlab 23.2.0.2428915 (R2023b) Update 4 for MACA64.  Copyright 2023 John J. Lee.
    
    properties
        testObj
    end
    
    methods (Test)
        function test_afun(this)
            import mlaif.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end
        function test_build_aif(this)
            %% Assess 3-5 sec frames vs 1 sec frames.
            %  Assess 5%, 10%, 50%, and 100% of centerline samples.
            %  Assess CO, OO, HO, FDG.

            
        end
    end
    
    methods (TestClassSetup)
        function setupMipIdif(this)
            import mlaif.*
            this.testObj_ = MipIdif();
        end
    end
    
    methods (TestMethodSetup)
        function setupMipIdifTest(this)
            this.testObj = this.testObj_;
            this.addTeardown(@this.cleanTestMethod)
        end
    end
    
    properties (Access = private)
        testObj_
    end
    
    methods (Access = private)
        function cleanTestMethod(this)
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
