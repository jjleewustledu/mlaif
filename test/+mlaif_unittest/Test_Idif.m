classdef Test_Idif < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 26-Sep-2023 23:19:23 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlaif/test/+mlaif_unittest.
    %  Developed on Matlab 9.14.0.2337262 (R2023a) Update 5 for MACI64.  Copyright 2023 John J. Lee.
    
    properties
        f
        f1
        t
        testObj
    end
    
    methods (Test)
        function test_afun(this)
            import mlaif.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end
        function test_moving_average(this)
            import mlaif.Idif.moving_average;

            figure
            plot(this.f1)
            hold on;
            g1 = moving_average(this.f1);
            this.verifyEqual(length(g1), 110)
            plot(g1, 'o:')
        end
        function test_moving_aver_oper(this)
            import mlaif.Idif.moving_aver_oper;
            import mlaif.Idif.moving_average;

            figure
            plot(this.f1, LineWidth=2)
            hold on;
            L = moving_aver_oper();
            this.verifyEqual(size(L), [110, 119])
            g2 = asrow(L*ascol(this.f1));
            plot(g2, 'o', MarkerSize=10)

            g1 = moving_average(this.f1);
            plot(g1, LineWidth=2)
            this.verifyEqual(g2, g1, RelTol=1e-15)

            g3 = asrow(lsqnonneg(L, ascol(g2)));
            plot(g3, 'o', MarkerSize=10)
            legend({'f1', 'L*f1', 'moving_average', 'lsqnonneg(L,g2)'}, Interpreter="none", FontSize=12)
        end
    end
    
    methods (TestClassSetup)
        function setupIdif(this)
            this.t = 0:118;
            this.f = exp(-this.t/10);
            this.f1 = zeros(size(this.f));
            this.f1(6:end) = this.f(1:end-5);
        end
    end
    
    methods (TestMethodSetup)
        function setupIdifTest(this)
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
