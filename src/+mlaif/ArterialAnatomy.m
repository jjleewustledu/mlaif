classdef (Abstract) ArterialAnatomy < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %% line1
    %  line2
    %  
    %  Created 07-Mar-2022 21:57:33 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlaif/src/+mlaif.
    %  Developed on Matlab 9.11.0.1873467 (R2021b) Update 3 for MACI64.  Copyright 2022 John J. Lee.
    
    methods (Static)
        function this = createForT1w(varargin)
            %% CREATE mlaif.ECIC using this factory method.  ECIC uses MP-RAGE.
            %  Args:
            %      bids (mlpipeline.Bids)

            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'bids', [], @(x) isa(x, 'mlpipeline.IBids'));
            parse(ip, varargin{:});
            ipr = ip.Results;
            
            this = mlaif.ECIC('anatomy', ipr.bids.t1w_ic, varargin{:});
        end
        function this = createForTof(varargin)
            %% CREATE mlaif.ICIC using this factory method.  ECIC uses TOF.
            %  Args:
            %      bids (mlpipeline.Bids)

            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'bids', [], @(x) isa(x, 'mlpipeline.IBids'));
            parse(ip, varargin{:});
            ipr = ip.Results;
            
            this = mlaif.ICIC('anatomy', ipr.bids.tof_ic, varargin{:});
        end
        function this = createForTofOnT1w(varargin)
            %% CREATE mlaif.ICIC using this factory method.  ECIC uses TOF.
            %  Args:
            %      bids (mlpipeline.Bids)

            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'bids', [], @(x) isa(x, 'mlpipeline.IBids'));
            parse(ip, varargin{:});
            ipr = ip.Results;
            
            this = mlaif.ICIC('anatomy', ipr.bids.tof_on_t1w_ic, varargin{:});
        end

        function fn = niigz(obj)
            if isa(obj, 'mlfourd.ImagingContext2')
                fn = strcat(obj.fqfp, '.nii.gz');
                return
            end
            ic = mlfourd.ImagingContext2(obj);
            fn = strcat(ic.fqfp, '.nii.gz');
        end
        function fn = mat(obj)
            if isa(obj, 'mlfourd.ImagingContext2')
                fn = strcat(obj.fqfp, '.mat');
                return
            end
            ic = mlfourd.ImagingContext2(obj);
            fn = strcat(ic.fqfp, '.mat');
        end
        function fn = json(obj)
            if isa(obj, 'mlfourd.ImagingContext2')
                fn = strcat(obj.fqfp, '.json');
                return
            end
            ic = mlfourd.ImagingContext2(obj);
            fn = strcat(ic.fqfp, '.json');
        end
    end

    properties (Dependent)
        anatomy
        anatPath
        anatTag
        bids
        dx
        dy
        dz
        imdilate_scale_mm % mm internaly converted to voxels and used with imdilate_strel
        imdilate_strel % default is 'sphere'
        mmppix
        Nx
        Ny
        Nz
        wmparc
    end

    methods

        %% GET/SET

        function g = get.anatomy(this)
            g = copy(this.anatomy_);
        end
        function g = get.anatPath(this)
            g = this.bids_.anatPath;
        end
        function g = get.anatTag(this)
            g = this.anatTag_;
        end
        function g = get.bids(this)
            g = this.bids_;
        end
        function g = get.dx(this)
            g = this.anatomy.imagingFormat.mmppix(1);
        end
        function g = get.dy(this)
            g = this.anatomy.imagingFormat.mmppix(2);
        end
        function g = get.dz(this)
            g = this.anatomy.imagingFormat.mmppix(3);
        end
        function g = get.imdilate_scale_mm(this)
            g = this.imdilate_scale_mm_;
        end
        function     set.imdilate_scale_mm(this, s)
            assert(isscalar(s) || isempty(s));
            this.imdilate_scale_mm_ = s;
        end
        function g = get.imdilate_strel(this)
            g = this.imdilate_strel_;
        end
        function     set.imdilate_strel(this, s)
            assert(istext(s));
            this.imdilate_strel_ = s;
        end
        function g = get.mmppix(this)
            if ~isempty(this.mmppix_)
                g = this.mmppix_;
                return
            end
            this.mmppix_ = this.anatomy.imagingFormat.mmppix;
            g = this.mmppix_;
        end
        function g = get.Nx(this)
            g = size(this.anatomy, 1);
        end
        function g = get.Ny(this)
            g = size(this.anatomy, 2);
        end
        function g = get.Nz(this)
            g = size(this.anatomy, 3);
        end
        function g = get.wmparc(this)
            g = this.wmparc_;
        end

        %%

        function [ic,f] = anatomy_on_pet(this, pet)
            %  Args:
            %      pet (any)
            %  Returns:
            %      ic (mlfourd.ImagingContext2): anatomy on pet.
            %      f (mlfsl.Flirt): usable in downstream flirt operations.

            pet = mlfourd.ImagingContext2(pet);

            [~,f_] = this.pet_static_on_anatomy(pet);

            tracer_tag = lower(this.bids.obj2tracer(pet));
            fqfn = strcat(this.anatomy.fqfp, '_on_', tracer_tag, '.nii.gz');
            ic = mlfourd.ImagingContext2(fqfn);
            ic = mlaif.Fung2013.sourcedata2derivatives(ic);
            f = copy(f_);
            f.invertXfm();
            f.in = this.niigz(this.anatomy);
            f.out = this.niigz(ic);
            f.ref = this.niigz(pet);
            f.applyXfm()
            assert(isfile(ic.fqfn))
        end
        function [ic,f] = pet_static_on_anatomy(this, pet)
            %  Args:
            %      pet (any): dyn is time-averaged to static as needed.
            %  Returns:
            %      ic (mlfourd.ImagingContext2): pet on anatomy.
            %      f (mlfsl.Flirt): usable in downstream flirt operations.

            pet = mlfourd.ImagingContext2(pet);
            if ~contains(pet.fileprefix, '_avgt') && ~contains(pet.fileprefix, 'static')
                pet = pet.timeAveraged();
            end
            pet.relocateDerivatives();

            if ~isfile(pet.fqfn)
                save(pet);
            end
            fqfn = strcat(pet.fqfp, '_on_', this.anatTag, '.nii.gz');
            ic = mlfourd.ImagingContext2(fqfn);
            f = mlfsl.Flirt( ...
                'in', pet, ...
                'ref', this.anatomy, ...
                'out', ic, ...
                'omat', this.mat(ic), ...
                'bins', 256, ...
                'cost', 'mutualinfo', ...
                'dof', 6, ...
                'interp', 'trilinear', ...
                'noclobber', true);
            f.flirt();
            assert(isfile(ic.fqfn))
        end
        function ic = pet_dyn_on_anatomy(this, pet_dyn, f)
            %  Args:
            %      pet_dyn (any): R^{3+1}l
            %      f (mlfsl.Flirt): from pet_static_on_anatomy().
            %  Returns:
            %      ic (mlfourd.ImagingContext2): pet_dyn on anatomy.

            pet_dyn = mlfourd.ImagingContext2(pet_dyn);
            assert(ndims(pet_dyn) == 4);
            assert(isa(f, 'mlfsl.Flirt'));
            pet_dyn.relocateDerivatives();
            
            if ~isfile(pet_dyn.fqfn)
                save(pet_dyn);
            end
            fqfn = strcat(pet_dyn.fqfp, '_on_', this.anatTag, '.nii.gz');
            ic = mlfourd.ImagingContext2(fqfn);
            f.in = pet_dyn;
            f.out = ic;
            f.applyXfm();
            assert(isfile(ic.fqfn));
        end  
    end

    %% PROTECTED

    properties (Access = protected)
        anatomy_
        anatTag_
        bids_
        imdilate_scale_mm_
        imdilate_strel_
        mmppix_
        wmparc_
    end

    methods (Access = protected)
        function this = ArterialAnatomy(varargin)
            %% ARTERIALANATOMY 
            %  Args:
            %      anatomy (any): understandable by mlfourd.ImagingContext2, specifies anatomical imaging, 
            %                     e.g., t1w or tof.
            %      bids (mlpipeline.Bids)
            %      imdilate_scale_mm (scalar):  mm of imdilate to apply to centerline before sampling PET.
            %      imdilate_strel (text):  structural element for imdilate; default is 'sphere'.
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            addParameter(ip, 'anatomy', []);
            addParameter(ip, 'bids', [], @(x) isa(x, 'mlpipeline.IBids'));
            addParameter(ip, 'imdilate_scale_mm', 0.25, @isscalar);
            addParameter(ip, 'imdilate_strel', 'sphere', @istext);
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            this.bids_ = ipr.bids;
            this.anatomy_ = mlfourd.ImagingContext2(ipr.anatomy);
            this.anatomy_ = this.bids.prepare_derivatives(this.anatomy_);
            this.anatTag_ = 'anat';
            this.imdilate_scale_mm_ = ipr.imdilate_scale_mm;
            this.imdilate_strel_ = ipr.imdilate_strel;
            this.wmparc_ = mlsurfer.Wmparc.createCoregisteredFromBids(this.bids_, this.anatomy_);
        end

        function that = copyElement(this)
            %%  See also web(fullfile(docroot, 'matlab/ref/matlab.mixin.copyable-class.html'))
            
            that = copyElement@matlab.mixin.Copyable(this);
            that.anatomy_ = copy(this.anatomy_);
            that.bids_ = copy(this.bids_);
            that.wmparc_ = copy(this.wmparc_);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
