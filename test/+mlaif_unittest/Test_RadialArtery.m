classdef Test_RadialArtery < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 28-Apr-2023 14:34:29 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlaif/test/+mlaif_unittest.
    %  Developed on Matlab 9.14.0.2239454 (R2023a) Update 1 for MACI64.  Copyright 2023 John J. Lee.
    
    properties
        aif % table with vars:  times, activityDensity
        testObj % mlaif.RadialArtery
        units % struct
        workdir = "/Users/jjlee/Library/CloudStorage/Box-Box/Matt Brier"
    end
    
    methods (Test)
        function test_afun(this)
            import mlaif.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end
        function test_ctor(this)
            disp(this.testObj)
        end
        function test_solve(this)
            obj = this.testObj.solve();
            plot(obj)
        end
        function test_writetable(this)
            obj = this.testObj.solve();
            writetable(obj, fullfile(this.workdir, "test_writetable.csv"));
            T = readtable(fullfile(this.workdir, "test_writetable.csv"));
            %disp(T)

            figure;
            plot(this.aif.times, this.aif.activityDensityDC, 'o', T.times, T.activityDensityDC)
        end
        function test_call(this)
            cd(this.workdir)
            for g = glob('JJL_1*.csv')'                
                aif_ori = readtable(g{1});
                aif_ = aif_ori;
                aif_.Properties.VariableNames = {'times', 'activityDensityDC'};
                aif_.times = floor(60*aif_.times);
                factor = 2.^(-aif_.times/6582);
                aif_.activityDensity = 37000 * aif_.activityDensityDC .* factor;
                obj = mlaif.RadialArtery( ...
                    'tracer', 'FDG', ...
                    'kernel', 1, ...
                    'model_kind', '3bolus', ...
                    'Measurement', aif_);
                obj = obj.solve();     
                the_csv = fullfile(this.workdir, strrep(g{1}, 'JJL_', 'JJL_mcmc_'));
                writetable(obj, the_csv);
                
                h = figure;
                T = readtable(the_csv);
                plot(aif_ori.Time_m_, aif_.activityDensityDC, 'o', T.times, T.activityDensityDC)
                title(the_csv, Interpreter='none');
                ylabel('decay-corr. activity (\muC/mL)', FontSize=18)
                xlabel('time (min)', FontSize=18)
                savemyfig(h, strrep(the_csv, '.csv', ''))
            end
        end
    end
    
    methods (TestClassSetup)
        function setupRadialArtery(this)
            import mlaif.*  
            this.aif = readtable(fullfile(this.workdir, "JJL_1179206_a.csv"));
            this.aif.Properties.VariableNames = {'times', 'activityDensityDC'};
            %this.aif.times = floor(60*this.aif.times);
            %factor = 2.^(-this.aif.times/6582);
            %this.aif.activityDensity = 37000 * this.aif.activityDensityDC .* factor;
            this.units.times = 'min';
            this.units.activityDensityDC = 'uCi/mL';

            this.testObj_ = RadialArtery( ...
                this.aif, ...
                this.units, ...
                tracer='FDG', ...
                model_kind='3bolus', ...
                kernel=1);
        end
    end
    
    methods (TestMethodSetup)
        function setupRadialArteryTest(this)
            if ishandle(this.testObj_)
                this.testObj = copy(this.testObj_);
            else
                this.testObj = this.testObj_;
            end
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
