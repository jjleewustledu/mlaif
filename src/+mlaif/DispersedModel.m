classdef DispersedModel < handle & mlaif.ArteryLee2021Model
    %% line1
    %  line2
    %  
    %  Created 03-Dec-2023 23:37:49 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlaif/src/+mlaif.
    %  Developed on Matlab 23.2.0.2428915 (R2023b) Update 4 for MACA64.  Copyright 2023 John J. Lee.
    
    methods
        function this = DispersedModel(varargin)
        end
    end

    methods (Static)
        function qs = sampled(ks, Data, ~, times_sampled)
            %% @return the Bayesian estimate of the measured AIF, including baseline, scaled to unity.
            
            N = Data.N;
            kernel = Data.kernel;
            baseline_frac = ks(9);
            scale_frac = 1 - baseline_frac;
            
            if kernel == 1
                W = 10;
                Data_ = Data;
                Data_.N = N+W-1;
                qs = mlaif.ArteryLee2021Model.solution(ks, Data_);
                qs = mlaif.ArteryLee2021Model.move_window(qs, W=W);
                qs = qs/max(qs); % \in [0 1]
                qs = scale_frac*qs + baseline_frac; % \in [0 1]
                return
            end

            qs = mlaif.ArteryLee2021Model.solution(ks, Data);
            if kernel ~= 1
                qs = conv(qs, kernel);
            end
            qs = qs(1:N);
            qs = qs/max(qs); % \in [0 1] 
            qs = scale_frac*qs + baseline_frac; % \in [0 1]   
                
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
