classdef ArterialCenterline < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %% ARTERIALCENTERLINE
    %  
    %  Created 06-Mar-2022 23:30:46 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlaif/src/+mlaif.
    %  Developed on Matlab 9.11.0.1873467 (R2021b) Update 3 for MACI64.  Copyright 2022 John J. Lee.

    properties (Dependent)
        anatomy % mlfourd.ImagingContext2(T1w|tof) aids visualization of centerline
        arterialAnatClass % text for logging
        cache_file % .mat
        coord1
        coord2
        pc_threshp
        product % interpolates pointCloud_for_centerline() to segmentation voxel barycenters
        segmentation
        suffix
        use_cache

        %% for B-splines in mlvg.Hunyadi2021
        
        k % order of b-splines ~ num. of inflections + 2
        t % represents knot positions in non-decreasing order with degeneracies representing higher derivatives
        M % grid coords, 3 x numel(grid)
        C % b-spline curve on neurologic grid, 3 x N_centerline_samples; see also bspline_deboor
        P % b-spline control points on neurologic grid as matrix ~ 3 x (len(t) - k); see alo bspine_estimate
        N_centerline_samples  % num. samples for bspline_deboor()
    end

    methods

        %% GET

        function g = get.anatomy(this)
            g = copy(this.fung2013_.anatomy);
        end
        function g = get.arterialAnatClass(this)
            g = this.fung2013_.arterialAnatClass;
        end
        function g = get.cache_file(this)
            g = fullfile(this.anatomy.filepath, ...
                strcat(clientname(true, 2), this.suffix, '.mat'));
        end
        function g = get.coord1(this)
            g = this.fung2013_.coord1;
        end
        function g = get.coord2(this)
            g = this.fung2013_.coord2;
        end
        function g = get.pc_threshp(this)
            g = this.fung2013_.pc_threshp;
        end
        function g = get.product(this)
            g = copy(this.product_);
        end
        function g = get.segmentation(this)
            g = this.segmentation_;
        end
        function g = get.suffix(this)
            c1 = this.fung2013_.str_coord1;
            c2 = this.fung2013_.str_coord2;
            g = sprintf('_%s-%sto%s', this.arterialAnatClass, c1, c2);
        end
        function g = get.use_cache(this)
            g = this.use_cache_;
        end

        function g = get.k(this)
            g = this.fung2013_.k;
        end
        function g = get.t(this)
            g = this.fung2013_.t;
        end
        function g = get.M(this)
            img = this.segmentation.imagingFormat.img;
            img = flip(img, 1);
            [I,J,K] = ind2sub(size(img), find(img));
            g(1,:) = I';
            g(2,:) = J';
            g(3,:) = K';
        end
        function g = get.C(this)
            g = this.C_;
        end
        function g = get.P(this)
            g = this.P_;
        end
        function g = get.N_centerline_samples(this)
            g = this.fung2013_.N_centerline_samples;
        end

        %%

        function this = ArterialCenterline(varargin)
            %% ARTERIALCENTERLINE 
            %  Args:
            %      fung2013 (mlaif.Fung2013): director of Fung2013 implementation, reference
            %      segmentation (mlfourd.ImagingContext2): product of mlaif.ArterialSegmentation.
            %      use_cache (logical): default is true.
            
            ip = inputParser;
            addParameter(ip, "fung2013", [], @(x) isa(x, 'mlaif.Fung2013'));
            addParameter(ip, "segmentation", [], @(x) isa(x, 'mlfourd.ImagingContext2') || isempty(x));
            addParameter(ip, "use_cache", true, @islogical);
            parse(ip, varargin{:});
            ipr = ip.Results;
            this.fung2013_ = ipr.fung2013;
            this.segmentation_ = ipr.segmentation;
            this.segmentation_.relocateDerivatives();
            this.use_cache_ = ipr.use_cache;

            this.hunyadi_ = mlvg.Hunyadi2021(); % adds package bspline to Matlab's path

            if this.use_cache_ && isfile(this.cache_file)
                this.load();
            end
        end

        function this = build_centerline(this)
            %% Builds a centerline using mlvg.Hunyadi2021.
            %  Returns:
            %      this.C: b-spline curve, 3 x N_centerline_samples; see also bspline_deboor.
            %      this.P: b-spline control points as matrix; see alo bspine_estimate.
            %      this.product: interpolates pointCloud_for_centerline() to segmentation voxel barycenters.

            if ~isempty(this.product_)
                return
            end 

            % de Boor's method
            this.P_ = bspline_estimate(this.k, this.t, this.M); % double
            this.C_ = bspline_deboor(this.k, this.t, this.P, this.N_centerline_samples); % double
            C__ = this.C_;
            save(strcat('debug_', clientname(true, 2), '.mat'), 'C__');

            % assemble product; save cache_file
            this.product_ = this.segmentation.zeros;
            this.product_.fileprefix = strcat(clientname(true, 2), this.suffix);
            this.product_.setPointCloud(this.pointCloud_for_centerline());
            save(this);
            plot3(this);
            pcshow(this);
        end
        function c = complement_coord(this, c, idx)
            c = this.fung2013_.complement_coord(c, idx);
        end
        function this = load(this)
            cf = load(this.cache_file);
            this.C_ = cf.this.C_;
            this.hunyadi_ = cf.this.hunyadi_;
            this.P_ = cf.this.P_;
            this.product_ = cf.this.product_;
            this.segmentation_ = cf.this.segmentation_;
        end
        function ax = pcshow(this, varargin)
            h = figure('Position', [1 1 1200 1200]);
            hold on; 
            pcshow(this.anatomy.pointCloud('threshp', this.pc_threshp));
            pc = this.pointCloud_for_anatomy();
            ax = pcshow(pc.Location, 'm', 'MarkerSize', 3);
            hold off;
            view(-210, 15) % tuned for orientstd
            daspect([1 1 1]);

            %set(h, 'InvertHardCopy', 'off');
            %set(h, 'Color', [0 0 0]);
            %set(gca, 'Color', [0 0 0]);
            xlabel('x (mm)', 'FontSize', 14)
            ylabel('y (mm)', 'FontSize', 14)
            zlabel('z (mm)', 'FontSize', 14)
            set(gca, 'xcolor', [1 1 1])
            set(gca, 'ycolor', [1 1 1])
            set(gca, 'zcolor', [1 1 1])

            title(strcat(clientname(true, 2), this.suffix), 'Interpreter', 'none', 'FontSize', 14);

            this.savefig(h);
        end
        function p = plot3(this, varargin)
            h = figure;
            hold on;
            p(1) = plot3(this.M(1,:), this.M(2,:), this.M(3,:), 'b.',  'MarkerSize', 1, varargin{:}); % segmentation
            p(2) = plot3(this.P(1,:), this.P(2,:), this.P(3,:), 'c--', 'MarkerSize', 4, varargin{:}); % B-spline control points
            p(3) = plot3(this.C(1,:), this.C(2,:), this.C(3,:), 'm',   'MarkerSize', 4, varargin{:}); % B-spline curve
            hold off;
            view(3);
            daspect([1 1 1]);
            axis tight;

            %set(h, 'InvertHardCopy', 'off');
            set(h, 'Color', [0 0 0]);
            set(gca, 'Color', [0 0 0]);
            xlabel('i (voxel)', 'FontSize', 14)
            ylabel('j (voxel)', 'FontSize', 14)
            zlabel('k (voxel)', 'FontSize', 14)
            set(gca, 'xcolor', [1 1 1])
            set(gca, 'ycolor', [1 1 1])
            set(gca, 'zcolor', [1 1 1])
       
            legend('segmentation', 'control points', 'b-spline curve', ...
                'Location', 'Best', 'EdgeColor', 'w', 'TextColor', 'w', 'FontSize', 11);
            title(strcat(clientname(true, 2), this.suffix), 'Interpreter', 'none', 'FontSize', 14);

            this.savefig(h);
        end
        function pc = pointCloud_for_anatomy(this, varargin)
            mmppix = this.anatomy.imagingFormat.mmppix;
            loc = this.C';

            coords = [this.coord1; this.coord2]; 
            originof = min(coords);
            sizeof = max(coords) - min(coords);
            originof = this.complement_coord(originof, 1);
            originof(1) = originof(1) - sizeof(1) - 1; % offset in neurologic orient

            loc(:,1) = (originof(1) + loc(:,1))*mmppix(1);
            loc(:,2) = (originof(2) + loc(:,2))*mmppix(2);
            loc(:,3) = (originof(3) + loc(:,3))*mmppix(3); % in xyz
            pc = pointCloud(loc, varargin{:}); 
        end
        function pc = pointCloud_for_centerline(this, varargin)
            mmppix = this.segmentation.imagingFormat.mmppix;
            loc = this.C';
            loc(:,1) = loc(:,1)*mmppix(1);
            loc(:,2) = loc(:,2)*mmppix(2);
            loc(:,3) = loc(:,3)*mmppix(3); % in xyz
            pc = pointCloud(loc, varargin{:}); 
        end
        function save(this)
            save(this.cache_file, 'this');
        end
        function savefig(this, h)
            fqfp = fullfile(fileparts(this.cache_file), clientname(true, 3));
            saveas(h, sprintf('%s%s.fig', fqfp, this.suffix));
            saveas(h, sprintf('%s%s.png', fqfp, this.suffix));
            close(h);
        end
    end
    
    %% PROTECTED

    properties (Access = protected)
        C_
        fung2013_
        hunyadi_
        P_
        product_
        segmentation_
        use_cache_
    end

    methods (Access = protected)
        function that = copyElement(this)
            %%  See also web(fullfile(docroot, 'matlab/ref/matlab.mixin.copyable-class.html'))
            
            that = copyElement@matlab.mixin.Copyable(this);
            that.hunyadi_ = copy(this.hunyadi_);
            that.product_ = copy(this.product_);
            that.segmentation_ = copy(this.segmentation_);
        end
    end

    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
