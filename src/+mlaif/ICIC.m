classdef ICIC < handle & mlaif.ArterialAnatomy
    %% INTRACRANIALINTERNALCAROTID
    %  
    %  Created 04-Mar-2022 14:07:51 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlaif/src/+mlaif.
    %  Developed on Matlab 9.11.0.1873467 (R2021b) Update 3 for MACI64.  Copyright 2022 John J. Lee.
    
    properties (Dependent)
        exclusion_blur   % numerical, depends on ECIC | ICIC
        exclusion_init   % mlfourd.ImagingContext2
        exclusion_list   % cell
        exclusion_thresh % numerical
        inclusion_blur   % numerical, depends on ECIC | ICIC
        inclusion_init   % mlfourd.ImagingContext2
        inclusion_thresh % numerical
    end

    methods

        %% GET

        function g = get.exclusion_blur(this)
            g = 0;
        end
        function g = get.exclusion_init(this)
            wmp = this.wmparc;
            g = copy(wmp.select_all());

            g = copy(wmp.select_roi('brainstem+'));
            mnp = round(10./this.wmparc.product.imagingFormat.mmppix); % 10 mm dilation
            g = g.imdilate_bin(strel("cuboid", mnp));            
            g = g | copy(wmp.select_all());
        end
        function g = get.exclusion_list(this)
            g = {};
        end
        function g = get.exclusion_thresh(this)
            g = [];
        end
        function g = get.inclusion_blur(this)
            g = 0;
        end
        function g = get.inclusion_init(this)
            if ~isempty(this.inclusion_init_)
                g = copy(this.inclusion_init_);
                return
            end
            g = this.anatomy.blurred(0.25);
            g = g.threshp(33);
            mnp = round(1./this.wmparc.product.imagingFormat.mmppix); % 1 mm dilation
            this.inclusion_init_ = g.imdilate_bin(strel("cuboid", mnp));
            g = copy(this.inclusion_init_);
        end
        function g = get.inclusion_thresh(this)
            g = [];
        end

        %%

        function this = ICIC(varargin)
            this = this@mlaif.ArterialAnatomy(varargin{:});
            this.anatTag_ = 'tof';
        end
    end    

    %% protected

    properties (Access = protected)
        inclusion_init_
    end

    methods (Access = protected)
        function that = copyElement(this)
            %%  See also web(fullfile(docroot, 'matlab/ref/matlab.mixin.copyable-class.html'))
            
            that = copyElement@matlab.mixin.Copyable(this);
            that.inclusion_init_ = copy(this.inclusion_init_);
        end
    end

    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
