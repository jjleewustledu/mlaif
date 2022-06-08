classdef (Abstract) AbstractFung2013 < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
	%% ABSTRACTFUNG2013 provides abstractions and reusable implementations of  
    %  Edward K Fung and Richard E Carson.  Cerebral blood flow with [15O]water PET studies using 
    %  an image-derived input function and MR-defined carotid centerlines.  
    %  Phys. Med. Biol. 58 (2013) 1903–1923.  doi:10.1088/0031-9155/58/6/1903
    %  
    %  Created 03-Mar-2022 19:27:58 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlaif/src/+mlaif.
    %  Developed on Matlab 9.11.0.1873467 (R2021b) Update 3 for MACI64.  Copyright 2022 John J. Lee.
    
	properties (Abstract)
        N_centerline_samples % for representation of b-splines
    end

    methods (Abstract)
        buildCenterlines(this)
        buildCorners(this)
        buildSegmentation(this)
        call(this)
    end

    methods (Static)
        function this = create(varargin)
            ip = inputParser;
            addParameter(ip, 'bids', [], @(x) isa(x, 'mlpipeline.Bids'));
            addParameter(ip, 'anatomy', [], @(x) isa(x, 'mlfourd.ImagingContext2'));
            parse(ip, varargin{:});
            this = [];
        end
    end

    properties 
        bbBuffer % extra voxels padded to coords to create convex bounding box for segmentation & centerline
        bbRange % coord ranges {x1:xN, y1:yN, z1:zN} for bounding box
        centerlines_ics % L, R in cell
        centerlines_pcs % L, R in cell
        contractBias % used by activecontour
        coords % 4 coord points along carotid centerlines at corners
        coords_b1_ic % coord points, blurred by 1 voxel fwhm, as ImagingContext2
        idifmask_ic % contains centerline modulated by innerRadius & outerRadius, as ImagingContext
        innerRadius % voxels; for sampling annulus
        iterations
        outerRadius % voxels; Fung reported best results with radius ~ 2.5, but implementation with strel requires integer num. of voxels.
        plotclose % close plots after saving
        plotdebug % show debugging plots
        ploton % show final results
        plotqc % show more plots for QA
        segmentation_blur
        segmentation_only
        segmentation_ic % contains solid 3D volumes for carotids
        segmentationThresh
        smoothFactor
        taus % containers.Map
        threshqc
        times % containers.Map
        timesMid % containers.Map 

        %% for B-splines in mlvg.Hunyadi2021

        k = 4
        t = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
        U % # samples for bspline_deboor()
        Cs % curves in cell, 3 x this.U
        Ps % control points in cell, 3 x M
    end

    properties (Dependent)
        anatomy
        anatomy_mask
        anatPath
        destinationPath
        petPath
        
        bids
    end

    methods

        %% GET, SET

        function g = get.anatomy(this)
            g = this.anatomy_;
        end
        function g = get.anatomy_mask(this)
            g = this.anatomy_mask_;
        end
        function g = get.anatPath(this)
            g = this.bids_.anatPath;
        end
        function g = get.bids(this)
            g = this.bids_;
        end
        function g = get.destinationPath(this)
            g = this.bids_.destinationPath;
        end
        function g = get.petPath(this)
            g = this.bids_.petPath;
        end

        %%

        function [pc,C,P] = buildCenterline(this, img, tag)
            %% Builds a centerline using mlvg.Hunyadi2021.
            %  If not previously saved,
            %  as requested by plotqc and plotdebug, plots centerline with thresholded anatomy and,
            %  as requested by plotclose, closes figures.
            %  
            %  Args:
            %      img (numeric): is the data upon which a centerline is built.
            %      tag (text): tags filenames.
            %  Returns:
            %      pc: is the pointCloud representation of the centerline.
            %      C: are points of the B-spline curve.
            %      P: is the matrix of B-spline control points.

            assert(istext(tag))
            idx = find(img);
            [X,Y,Z] = ind2sub(size(img), idx);             
            M(1,:) = X'; % M are ints cast as double
            M(2,:) = Y';
            M(3,:) = Z';
            this.U = this.N_centerline_samples;             
            P = bspline_estimate(this.k, this.t, M); % double
            C = bspline_deboor(this.k, this.t, P, this.U); % double, ~2x oversampling for Z
            pc = pointCloud(C');
            
            class_str = strrep(class(this), 'mlaif.', '');
            fp = fullfile(this.destinationPath, ...
                sprintf('%s_%s_centerline_in_%s', class_str, tag, this.anatomy.fileprefix));
            if this.plotqc
                h = figure;
                pcshow(pointCloud(this.anatomy, 'thresh', this.threshqc*dipmax(this.anatomy)))
                hold on; pcshow(pc.Location, '*m', 'MarkerSize', 12); hold off;
                saveas(h, [fp '.fig'])
                set(h, 'InvertHardCopy', 'off');
                set(h,'Color',[0 0 0]); % RGB values [0 0 0] indicates black color
                saveas(h, [fp '.png'])
                if this.plotclose
                    close(h)
                end
            end
            fp1 = fullfile(this.destinationPath, ...
                sprintf('%s_%s_centerline_in_segmentation', class_str, tag));
            if this.plotdebug
                h1 = figure;
                hold all;
                plot3(M(1,:), M(2,:), M(3,:), 'g.');
                plot3(P(1,:), P(2,:), P(3,:), 'b');
                plot3(C(1,:), C(2,:), C(3,:), 'r');
                legend('segmentation', 'control points', 'curve', ...
                    'Location', 'Best');
                hold off;
                saveas(h1, [fp1 '.fig'])
                saveas(h1, [fp1 '.png'])
                if this.plotclose
                    close(h1)
                end
            end
        end
        function h = plotIdif(this, tbl_idif)
            %% As requested by ploton, plots then saves all IDIFs in the subject collection.  Clobbers previously saved.
            %  As requested by plotclose, closes figures.

            if this.ploton
                h = figure;
                hold on
                tracer_ = tbl_idif.tracer;
                for irow = 1:size(tbl_idif,1)
                    timesMid_ = this.timesMid(tracer_{irow});
                    IDIF_ = tbl_idif.IDIF{irow};
                    N = min(length(timesMid_), length(IDIF_));
                    switch tracer_{irow}
                        case {'OC' 'CO'}
                            linestyle = '-.';
                        case 'OO'
                            linestyle = '-';
                        case 'HO'
                            linestyle = '--';
                        otherwise
                            linestyle = ':';
                    end
                    plot(timesMid_(1:N), IDIF_(1:N), linestyle)
                end
                xlim([0 350])
                xlabel('time (s)')
                ylabel('activity density (Bq/mL)')
                title('Image-derived Input Functions')
                legend(tracer_')
                hold off
                [~,fqfp] = fileparts(tbl_idif.Properties.Description);
                saveas(h, [fqfp '.fig'])
                saveas(h, [fqfp '.png'])
                if this.plotclose
                    close(h)
                end
            end
        end
        function h = plotSegmentation(this, ac, varargin)
            %% As requested by plotqc, plots then saves segmentations by activecontour.
            %  As requested by plotclose, closes figure.
            %  @param required activecontour result.
            %  @param iterations is integer.
            %  @param smoothFactor is scalar.

            ip = inputParser;
            addRequired(ip, 'ac', @islogical)
            addOptional(ip, 'iterations', this.iterations, @isscalar)
            addOptional(ip, 'smoothFactor', this.smoothFactor, @isscalar)
            parse(ip, ac, varargin{:})
            ipr = ip.Results;
            
            fp = fullfile(this.destinationPath, [this.anatomy.fileprefix '_snakes']);
            if this.plotqc
                h = figure;
                mmppix = this.anatomy.imagingFormat.mmppix;
                L = (size(ipr.ac) - 1) .* mmppix;
                [X,Y,Z] = meshgrid(0:mmppix(1):L(1), 0:mmppix(2):L(2), 0:mmppix(3):L(3));
                X = permute(X, [2 1 3]);
                Y = permute(Y, [2 1 3]);
                Z = permute(Z, [2 1 3]);
                p = patch(isosurface(X, Y, Z, double(ipr.ac)));
                p.FaceColor = 'red';
                p.EdgeColor = 'none';
                daspect([1 1 1])
                camlight;
                lighting phong
                title(sprintf('iterations %i, contractBias %g, smoothFactor %g, segmentationThresh %g', ...
                    ipr.iterations, this.contractBias, ipr.smoothFactor, this.segmentationThresh))
                saveas(h, [fp '.fig'])
                saveas(h, [fp '.png'])
                if this.plotclose
                    close(h)
                end
            end
        end
        function writetable(this, t, activity, dynfp)
            len = min(length(t), length(activity));
            t = ascol(t(1:len));
            activity = ascol(activity(1:len));
            tbl = table(t, activity);
            tbl.Properties.Description = [class(this) '_' this.subjectFolder];
            tbl.Properties.VariableUnits = {'s', 'Bq/mL'};

            fqfn = fullfile(this.destinationPath, sprintf('%s_idif.csv', dynfp));
            writetable(tbl, fqfn)
        end
    end

    %% PROTECTED

    properties (Access = protected)
        anatomy_
        anatomy_mask_
        bids_
        hunyadi_
        sessionData_
    end
    methods (Access = protected)
        function this = AbstractFung2013(varargin)
            %% ABSTRACTFUNG2013 
            %  Args:
            %      arg1 (its_class): Description of arg1.
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'bids', [], @(x) isa(x, 'mlpipeline.Bids'))
            addParameter(ip, 'contractBias', 0.02, @(x) isscalar(x) || isa(x, 'containers.Map'))
            addParameter(ip, 'corners', [], @(x) ismatrix(x) || isa(x, 'containers.Map'))
            addParameter(ip, 'innerRadius', 0, @isnumeric)
            addParameter(ip, 'iterations', 70, @(x) isscalar(x) || isa(x, 'containers.Map'))
            addParameter(ip, 'outerRadius', 2, @isnumeric)
            addParameter(ip, 'ploton', true, @islogical)
            addParameter(ip, 'plotqc', true, @islogical)
            addParameter(ip, 'plotdebug', false, @islogical)
            addParameter(ip, 'plotclose', true, @islogical)
            addParameter(ip, 'segmentationBlur', 0, @(x) isscalar(x) || isa(x, 'containers.Map'))
            addParameter(ip, 'segmentationOnly', false, @islogical)
            addParameter(ip, 'segmentationThresh', 190, @(x) isscalar(x) || isa(x, 'containers.Map'))
            addParameter(ip, 'sessionData', [], @(x) isa(x, 'mlpipeline.ISessionData') || isempty(x))
            addParameter(ip, 'smoothFactor', 0, @isscalar)
            addParameter(ip, 'threshqc', 0.75, @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            if isa(ipr.contractBias, 'containers.Map')
                ipr.contractBias = ipr.contractBias(this.subjectFolder);
            end 
            if isa(ipr.corners, 'containers.Map')
                ipr.corners = ipr.corners(this.subjectFolder);
            end
            if isa(ipr.iterations, 'containers.Map')
                ipr.iterations = ipr.iterations(this.subjectFolder);
            end
            if isa(ipr.segmentationBlur, 'containers.Map')
                ipr.segmentationBlur = ipr.segmentationBlur(this.subjectFolder);
            end
            if isa(ipr.segmentationThresh, 'containers.Map')
                ipr.segmentationThresh = ipr.segmentationThresh(this.subjectFolder);
            end
            this.bids_ = ipr.bids;
            this.contractBias = ipr.contractBias;
            this.coords = ipr.corners;
            this.innerRadius = ipr.innerRadius;
            this.iterations = ipr.iterations;
            this.outerRadius = ipr.outerRadius;
            this.ploton = ipr.ploton;
            this.plotqc = ipr.plotqc;
            this.plotdebug = ipr.plotdebug;
            this.plotclose = ipr.plotclose;
            this.segmentation_blur = ipr.segmentationBlur;
            this.segmentation_only = ipr.segmentationOnly;
            this.segmentationThresh = ipr.segmentationThresh;
            this.sessionData_ = ipr.sessionData;
            this.smoothFactor = ipr.smoothFactor;
            this.threshqc = ipr.threshqc;

            this.hunyadi_ = mlvg.Hunyadi2021();
        end

        function that = copyElement(this)
            %%  See also web(fullfile(docroot, 'matlab/ref/matlab.mixin.copyable-class.html'))
            
            that = copyElement@matlab.mixin.Copyable(this);
            that.bids_ = copy(this.bids_);
            that.hunyadi_ = copy(this.hunyadi_);
        end
        function decay_uncorrected = decay_uncorrected(this, idif)
            %  @param idif is an mlfourd.ImagingContext2 containing a double row.
            %  @returns decay_uncorrected, the IDIF as a double row.

            assert(isa(idif, 'mlfourd.ImagingContext2'))
            decay_corrected = idif.nifti.img;
            assert(isvector(decay_corrected))

            tracer = this.bids.obj2tracer(idif);
            taus_ = this.taus(tracer);
            N = min(length(decay_corrected), length(taus_));
            radio = mlpet.Radionuclides(tracer);
            decay_uncorrected = decay_corrected(1:N) ./ radio.decayCorrectionFactors('taus', taus_(1:N));
            decay_uncorrected = asrow(decay_uncorrected);
        end
        function box = ensureBoxInFieldOfView(this, box)
            %% removes any elements of box := {xrange yrange zrange}, with anisotropic ranges, 
            %  that lie outside of the field of view of this.anatomy.

            assert(iscell(box), 'mlaif:ValueError', ...
                'AbstractFung2013.ensureBoxInFieldOfView: class(box)->%s', class(box))
            size_ = size(this.anatomy);
            for m = 1:length(box)
                bm = box{m};
                box{m} = bm(1 <= bm & bm <= size_(m));
            end
        end
        function [X,Y,Z] = ensureSubInFieldOfView(this, X, Y, Z)
            %% removes any subscripts X, Y, and Z, which are equally sized, that lie outside of 
            %  the field of view of this.anatomy.

            assert(isvector(X))
            assert(isvector(Y))
            assert(isvector(Z))
            assert(length(X) == length(Y) && length(Y) == length(Z))

            size_ = size(this.anatomy);
            toss_ =          X < 1 | size_(1) < X;
            toss_ = toss_ | (Y < 1 | size_(2) < Y);
            toss_ = toss_ | (Z < 1 | size_(3) < Z);

            X = X(~toss_);
            Y = Y(~toss_);
            Z = Z(~toss_);
        end
        function ic = pointCloudToIC(this, pc, varargin)
            ip = inputParser;
            addRequired(ip, 'pc', @(x) isa(x, 'pointCloud'))
            addOptional(ip, 'fileprefix', 'pointCloudToIC', @ischar)
            parse(ip, pc, varargin{:})
            ipr = ip.Results;
            
            ifc = this.anatomy.nifti;
            ifc.fileprefix = ipr.fileprefix;
            X = round(pc.Location(:,1));
            Y = round(pc.Location(:,2));
            Z = round(pc.Location(:,3));
            [X,Y,Z] = this.ensureSubInFieldOfView(X, Y, Z);
            ind = sub2ind(size(ifc), X, Y, Z);
            img = zeros(size(this.anatomy));
            img(ind) = 1;
            ifc.img = img;
            ic = mlfourd.ImagingContext2(ifc);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
