classdef ArterialBoundingBox < handle & matlab.mixin.Heterogeneous & matlab.mixin.Copyable
    %% Implements a bounding box for an arterial structure useful for 
    %  ArterialSegmentation, ArterialCenterline.
    %  Usage: >> coord1 = [x1, y1, z1]; % e.g. in RPI corner of visible arterial structure
    %         >> coord2 = [x1, y2, z2]; % e.g. in LAS corner of structure
    %         >> abb = mlaif.ArterialBoundingBox(imgobj, coord1, coord2);
    %         >> e = abb.extract(); % bounding box containing structure of interest
    %         >> e = e * 3; % brightened e, but can by any grid-preserving transformation
    %         >> i = abb.insert(e); % recovers imgobj with brightening of e
    %
    %  Created 22-Mar-2022 23:24:40 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlaif/src/+mlaif.
    %  Developed on Matlab 9.12.0.1884302 (R2022a) for MACI64.  Copyright 2022 John J. Lee.
    
    properties (Dependent)
        imagingContext % full field of view
        mmppix % of imagingContext
        originof % of bounding box in voxel indices
        sizeof % of bounding box in voxel indices
    end

    methods

        %% GET

        function g = get.imagingContext(this)
            g = copy(this.imagingContext_);
        end
        function g = get.mmppix(this)
            g = this.imagingContext_.imagingFormat.mmppix;
        end
        function g = get.originof(this)
            g = this.originof_;
        end
        function g = get.sizeof(this)
            g = this.sizeof_;
        end

        %% 

        function this = enlarge_box_base(this, varargin)
            %% Enlarges box by scalar multiple in x-y plane
            %  Args:
            %      multiple (scalar): e.g., multiple->3 enlarges box of size [10 10 100] to [30 30 100].

            ip = inputParser;
            addOptional(ip, 'multiple', 3, @isscalar);
            parse(ip, varargin{:});
            ipr = ip.Results;

            D = round(ipr.multiple*this.sizeof_(1:2));
            Dshift = round(D/2 - this.sizeof_(1:2)/2);
            this.originof_(1:2) = this.originof_(1:2) - Dshift; 
            this.sizeof_(1:2) = D;
        end
        function e = extract(this)
            %  Returns:
            %      e (mlfourd.ImagingContext2): from this.imagingContext with size(e) ~ this.sizeof.

            this.imagingContext_.relocateToDerivativesFolder();

            o = this.originof_;
            s = this.sizeof_;
            if ~isfile(this.imagingContext_.fqfn)
                this.imagingContext_.save();
            end
            e = this.imagingContext_.zoomed(o(1), s(1), o(2), s(2), o(3), s(3));
        end
        function ic = insert(this, ins)
            %  Args:
            %      ins (any): understood by mlfourd.ImagingContext2; size(ins) must cover this.sizeof.
            %  Returns:
            %      ic (mlfourd.ImagingContext2): this.imagingContext with insertion at this.origin, this.size;
            %                                    ins is cropped to this.size as needed.

            ins = mlfourd.ImagingContext2(ins);
            ins.relocateToDerivativesFolder();

            if 3 == ndims(ins)
                ic = this.insert3(ins);
                return
            end
            if 4 == ndims(ins)
                ic = this.insert4(ins);
                return
            end
            error('mlaif:NotImplementedError', 'ArterialBoundingBox.insert(): ndims(ins)->%i', ndims(ins));
        end

        function this = ArterialBoundingBox(varargin)
            %% ARTERIALBOUNDINGBOX 
            %  Args:
            %      imgobj (required, any): understood by mlfourd.ImagingContext2.
            %      coord1 (required, numeric): voxel i,j,k for exterior of 1st corner.
            %      coord2 (required, numeric): voxel i,j,k for interior of 2nd corner, diagonal from 1st.
            
            ip = inputParser;
            addRequired(ip, "imgobj")
            addRequired(ip, "coord1", @isnumeric)
            addRequired(ip, "coord2", @isnumeric)
            parse(ip, varargin{:})
            ipr = ip.Results;
            
            this.imagingContext_ = mlfourd.ImagingContext2(ipr.imgobj);
            coords = [asrow(ipr.coord1); asrow(ipr.coord2)]; % 2 x 3
            this.originof_ = min(coords);
            this.sizeof_ = max(coords) - min(coords);
        end
    end

    %% PROTECTED

    properties (Access = protected)
        imagingContext_
        originof_
        sizeof_
    end

    methods (Access = protected)
        function that = copyElement(this)
            %%  See also web(fullfile(docroot, 'matlab/ref/matlab.mixin.copyable-class.html'))
            
            that = copyElement@matlab.mixin.Copyable(this);
            that.imagingContext_ = copy(this.imagingContext_);
        end
        function ic = insert3(this, ins)
            ins = mlfourd.ImagingContext2(ins);
            ins_img = ins.imagingFormat().img;

            o = this.originof_ + 1;
            d = o + this.sizeof_ - 1;
            ifc = this.imagingContext_.imagingFormat();
            if any(d > size(ifc))
                d = min([d; size(ifc)]); % trim ins to fit
                ifc.img(o(1):d(1), o(2):d(2), o(3):d(3)) = ins_img(1:d(1), 1:d(2), 1:d(3));
            else
                ifc.img(o(1):d(1), o(2):d(2), o(3):d(3)) = ins_img;
            end
            ic = mlfourd.ImagingContext2(ifc);
        end
        function ic = insert4(this, ins)
            ins = mlfourd.ImagingContext2(ins);
            ins_img = ins.imagingFormat().img;

            o = this.originof_ + 1;
            d = o + this.sizeof_ - 1;
            ifc = this.imagingContext_.imagingFormat();
            sz = size(ifc);
            if any(d > sz(1:3))
                d = min([d; sz(1:3)]); % trim ins to fit
                ifc.img(o(1):d(1), o(2):d(2), o(3):d(3), :) = ins_img(1:d(1), 1:d(2), 1:d(3), :);
            else
                ifc.img(o(1):d(1), o(2):d(2), o(3):d(3), :) = ins_img;
            end
            ic = mlfourd.ImagingContext2(ifc);
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
