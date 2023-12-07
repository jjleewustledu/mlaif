classdef Test_ArteryLee2021Model < matlab.unittest.TestCase
    %% line1
    %  line2
    %  
    %  Created 04-Dec-2023 00:14:58 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlaif/test/+mlaif_unittest.
    %  Developed on Matlab 23.2.0.2428915 (R2023b) Update 4 for MACA64.  Copyright 2023 John J. Lee.
    
    properties
        bcm
        testObj
    end
    
    methods (Test)
        function test_afun(this)
            import mlaif.*
            this.assumeEqual(1,1);
            this.verifyEqual(1,1);
            this.assertEqual(1,1);
        end

        %% BoxcarModel
        
        function test_BoxcarModel_build_simulated(this)
            ks = [0.5, 0.05, 1, 0, 0, 0.2, 0.05, 30, 0];
            plot(0:80, this.bcm.build_simulated(ks))
        end
        function test_BoxcarModel_build_solution(this)
            soln = this.bcm.build_solution();
            disp(soln)
            obj = this.bcm;
            figure; plot(obj.times_sampled, obj.artery.imagingFormat.img)
            %disp(obj.times_sampled)
            %disp(obj.artery_interpolated)
        end
        function test_BoxcarModel_sampled(this)
            ks = [0.5, 0.05, 1, 0, 0, 0.2, 0.05, 30, 0];
            Data = this.bcm.Data;
            times_sampled = this.bcm.times_sampled;
            qs = this.bcm.sampled(ks, Data, [], times_sampled);
            plot(times_sampled, qs, ':o');
        end
        function test_BoxcarModel_solution(this)
            ks = [0.5, 0.05, 1, 0, 0, 0.2, 0.05, 30, 0];
            Data = this.bcm.Data;
            times_sampled = this.bcm.times_sampled;
            qs = this.bcm.solution(ks, Data);
            plot(0:times_sampled(end), qs, ':o');
        end

        %% ArteryLee2021Model
        
        function test_build_simulated(this)
            ks = [0.5, 0.05, 1, 0, 0, 0.2, 0.05, 30, 0];
            plot(0:80, this.testObj.build_simulated(ks))
        end
        function test_build_solution(this)
            soln = this.testObj.build_solution();
            disp(soln)
            obj = this.testObj;
            figure; plot(obj.times_sampled, obj.artery.imagingFormat.img)
            %disp(obj.times_sampled)
            %disp(obj.artery_interpolated)
        end
        function test_sampled(this)
            ks = [0.5, 0.05, 1, 0, 0, 0.2, 0.05, 30, 0];
            Data = this.testObj.Data;
            times_sampled = this.testObj.times_sampled;
            qs = this.testObj.sampled(ks, Data, [], times_sampled);
            plot(times_sampled, qs, ':o');
        end
        function test_solution(this)
            ks = [0.5, 0.05, 1, 0, 0, 0.2, 0.05, 30, 0];
            Data = this.testObj.Data;
            times_sampled = this.testObj.times_sampled;
            qs = this.testObj.solution(ks, Data);
            plot(0:times_sampled(end), qs, ':o');
        end
    end
    
    methods (TestClassSetup)
        function setupArteryLee2021Model(this)
            ic = mlfourd.ImagingContext2( ...
                fullfile(getenv("HOME"), "MATLAB-Drive", "mlaif", "data", ...
                "sub-108293_ses-20210421150523_trc-oo_proc-delay0-BrainMoCo2-createNiftiMovingAvgFrames-ScannerKit-do-make-activity-density_timeAppend-4_pet_MipIdif_build_aif_cellTimingData.nii.gz"));
            this.testObj_ = mlaif.ArteryLee2021Model.create( ...
                artery=ic, tracer="OO", closeFigures=false, Nensemble=30);
            this.bcm_ = mlsiemens.BoxcarModel.create( ...
                artery=ic, tracer="OO", closeFigures=false, Nensemble=10);

            disp(ic.filepath)
        end
    end
    
    methods (TestMethodSetup)
        function setupArteryLee2021ModelTest(this)
            this.testObj = copy(this.testObj_);
            this.bcm = copy(this.bcm_);
            this.addTeardown(@this.cleanTestMethod)
        end
    end
    
    properties (Access = private)
        bcm_
        testObj_
    end
    
    methods (Access = private)
        function cleanTestMethod(this)
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
