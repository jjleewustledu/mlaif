classdef WholeBrainPET < mlaif.PETDynamicsAbstract
	%% WHOLEBRAINPET manages dynamic, masked PET .
    %  KLUDGE:  PET time samples hard-coded to be 0.5 Hz
	%
	%  $Revision$ 
	%  was created $Date$ 
	%  by $Author$, 
	%  last modified $LastChangedDate$ 
	%  and checked into repository $URL$, 
	%  developed on Matlab 8.3.0.532 (R2014a) 
	%  $Id$ 
	%  N.B. classdef (Sealed, Hidden, InferiorClasses = {?class1,?class2}, ConstructOnLoad) 

    properties (Constant)
        PET_NATIVE_TIMES = 2:2:120;
    end
    
    properties        
        nativeMeasurements
        nativeTimes
        tracerConcentrations
    end
    
    properties (Dependent)
        fileprefix
 	end 

    methods %% GET
        function fp = get.fileprefix(this)
            assert(~isempty(this.pet_));
            fp = this.pet_.fileprefix;
        end
    end
    
    methods (Static)  
        function this   = load(petfn, varargin)
            %% LOAD
            %  Usage:   this = WholeBrainPET.load(pet_filename[, parameter_name, parameter_value ...])

            this = mlaif.WholeBrainPET(petfn, varargin{:});
        end
        function          flipAllPET(pth)
            assert(lexist(pth, 'dir'));
            
            pnum = mlaif.WholeBrainPET.pNumber(pth);
            flipPETs(fullfile(pth, [    pnum '*.nii.gz']));            
            flipPETs(fullfile(pth, ['r' pnum '*.nii.gz']));
            flipPETs(fullfile(pth,         'cs*.nii.gz'));
            
            function flipPETs(str)
                dt   = mlsystem.DirTool(str);
                for f = 1:length(dt)
                   mlaif.WholeBrainPET.flipSinglePET(dt.fqfns{f});
                end  
            end
        end
        function          flipSinglePET(fname)
            nii     = mlfourd.NIfTI.load(fname);
            nii.img = flip4d(nii.img, 'y');
            nii.save;
        end
        function          maskAllDynamicPET
            import mlaif.*;
            
            dt = mlsystem.DirTools('mm*');
            for f = 1:dt.length
                try
                    pnum     = WholeBrainPET.pNumber(dt.fqdns{f});
                    maskfn   = fullfile(dt.fqdns{f}, 'fsl', 'bt1_default_mask_on_ho_meanvol_default.nii.gz');
                    petfn    = fullfile(dt.fqdns{f}, 'ECAT_EXACT', '962_4dfp', [pnum 'ho1.nii.gz']);
                    maskedfn = fullfile(dt.fqdns{f}, 'ECAT_EXACT', '962_4dfp', [pnum 'ho1_masked.nii.gz']);
                    WholeBrainPET.maskSingleDynamicPET(maskfn, petfn, maskedfn);
                    fprintf('wrote %s\n', maskedfn);
                catch ME 
                    handwarning(ME);
                    % fprintf('no dynamic PET found in %s\n', dt.fqdns{f});
                end
            end
        end
        function petfn1 = maskSingleDynamicPET(maskfn, petfn, petfn1)
            %% MASKDYNAMICPET reads a prepared mask and a dynamic PET study, then applies the mask and 
            %  saves the masked PET study; masked PET filename default is [fileprefix_original '_masked.nii.gz'].
            %  Usage:  pet_masked_filename = WholeBrainPET.maskDynamicPET(mask, pet_NIfTI_filename[, pet_masked_filename])
            
            import mlaif.*;
            if (~exist('petfn1', 'var'))
                petfn1 = [fileprefix(petfn) '_masked.nii.gz']; end

            import mlfourd.*;
            mask = NIfTI.load(maskfn);
            assert(logical(mask.dipmax <= 1));
            pet = NIfTI.load(petfn);
            WholeBrainPET.assertVolumesEqual(mask, pet);
            
            img = pet.img;
            for t = 1:size(pet,4)
                img(:,:,:,t) = img(:,:,:,t) .* mask.img;
            end
            pet.img = img;
            pet.saveas(petfn1);
        end        
        function          renameAllDynos
            import mlaif.*;
            
            pwd0 = pwd;
            dt = mlsystem.DirTools('mm*');
            
            for f = 1:dt.length
                targetFolder = fullfile(dt.fqdns{f}, 'ECAT_EXACT', '962_4dfp', '');
                if (lexist(targetFolder, 'dir'))
                    try
                        cd(targetFolder);
                        WholeBrainPET.renameSomeDynos;
                    catch ME %#ok<NASGU>
                        fprintf('PROBLEM with renameDynaFiles in folder %s\n', dt.fqdns{f});
                    end
                    cd(pwd0);
                end
            end
        end
        function          renameSomeDynos
            dt  = mlfourd.DirTool('*dyno*');
            for f = 1:dt.length
                baseFilename   = dyno2base(  dt.fns{f});
                summedFilename = base2summed(baseFilename);
                system(sprintf('mv %s %s',   baseFilename, summedFilename));
                system(sprintf('mv %s %s',   dt.fns{f},    baseFilename));                        
            end
            
            function b = dyno2base(d)                
                idx = strfind(d, '_dyno');
                b   = sprintf('%s.nii.gz', d(1:idx-1));
            end
            function s = base2summed(b)
                s = sprintf('%s_summed.nii.gz', b(1:end-7));
            end
        end
        function          assertVolumesEqual(nii, nii1)
            assert(all(niftiVolume(nii) == niftiVolume(nii1)));
            
            function vol = niftiVolume(nii)
                assert(isa(nii, 'mlfourd.NIfTI'));
                vol = size(nii);
                vol = vol(1:3);
            end
        end
        function fn     = tracerFilestem(fn)
            %% TRACERFILESTEM strips suffixes such as _fxtoy, .nii.gz
            %  e.g.,   WholeBrain.tracerFilename('p1234ho1_f9to99.nii.gz')
            %          p1234ho1
            
            fn = strtok(fn, '.');
            fn = strtok(fn, '_');
        end
        function pnum   = pNumber(str)
            %% PNUMBER returns the pnumber, p1234, from an arbitrary string using regexp            
            
            pnum  = str2pnum(str);
        end
    end
    
    %% PRIVATE
    
    properties (Access = 'private')
        pet_
    end
    
    
	methods (Access = 'private')
        
 		function this = WholeBrainPET(petfn, varargin)
			%% WHOLEBRAINPET 
 			%  Usage:   this = WholeBrainPET([parameter, parameter_values, ...])
            %                                ^ TimeInterpolants
            
            this = this@mlaif.PETDynamicsAbstract(varargin{:});
            
            this.pet_ = mlfourd.NIfTI.load(petfn);
            this.nativeMeasurements = squeeze(sum(sum(sum(this.pet_.img,1), 2), 3)); % integral dr pet(r)
            this.nativeTimes = this.PET_NATIVE_TIMES;            
            
            this.tracerConcentrations = this.clipFirstMeasurement(this.nativeMeasurements);
            if (~this.equivalent(this.timeInterpolants, this.nativeTimes))
                this.tracerConcentrations = this.spline; end
		end 
 	end 

	%  Created with NewClassStrategy by John J. Lee, after newfcn by Frank Gonzalez-Morphy 
end

