classdef SessionsCatalog  
	%% SESSIONSCATALOG   

	%  $Revision$ 
 	%  was created $Date$ 
 	%  by $Author$,  
 	%  last modified $LastChangedDate$ 
 	%  and checked into repository $URL$,  
 	%  developed on Matlab 8.4.0.150421 (R2014b) 
 	%  $Id$ 
 	 

    properties (Constant)
        FSL_FOLDER = 'fsl'
        PERFUSION_4DFP_FOLDER = 'perfusion_4dfp'
        PERFUSION_4DFP_LOG = 'perfusion_4dfp.log'
        CRV_FOLDER = 'ECAT_EXACT/pet'
        HO_FOLDER = 'ECAT_EXACT/962_4dfp'
    end
    
	properties 
 		 npPath
         hasEp2dDefaultMcf % session folder names, not path
    end
    
	methods 
        function fqfn = perfusion4dfpLog(this, f)
            fqfn = fullfile(this.npPath, this.hasEp2dDefaultMcf{f}, this.PERFUSION_4DFP_FOLDER, this.PERFUSION_4DFP_LOG);
        end
        function fqfn = crvFile(this, f)
            fqfn = fullfile(this.npPath, this.hasEp2dDefaultMcf{f}, this.CRV_FOLDER, [str2pnum(this.hasEp2dDefaultMcf{f}) 'ho1.crv']);
        end
        function fqfn = hoFile(this, f)            
            fqfn = fullfile(this.npPath, this.hasEp2dDefaultMcf{f}, this.HO_FOLDER, [str2pnum(this.hasEp2dDefaultMcf{f}) 'ho1.nii.gz']);
        end
        function fqfn = maskFile(this, f)            
            fqfn = fullfile(this.npPath, this.hasEp2dDefaultMcf{f}, this.FSL_FOLDER, 'bt1_default_mask_on_ho_meanvol_default_161616fwhh.nii.gz');
        end
        function fqfn = maskFile2(this, f)            
            fqfn = fullfile(this.npPath, this.hasEp2dDefaultMcf{f}, this.FSL_FOLDER, 'bt1_default_mask_on_ho_meanvol_default.nii.gz');
        end
        function fqfn = hoMasked(this, f)            
            fqfn = fullfile(this.npPath, this.hasEp2dDefaultMcf{f}, this.HO_FOLDER, [str2pnum(this.hasEp2dDefaultMcf{f}) 'ho1_masked.nii.gz']);
        end
        
        function m    = missingPerfusion4dfpLog(this)
            m = this.findSomethingMissing(@this.perfusion4dfpLog);
        end
        function m    = missingCrvFile(this) 
            m = this.findSomethingMissing(@this.crvFile);
        end
        function m    = missingHoFile(this)
            m = this.findSomethingMissing(@this.hoFile);
        end
        function m    = missingMaskFile(this)
            m = this.findSomethingMissing(@this.maskFile2);
        end
        function inc  = checkConsistencies(this)            
            inc = {};
            import mlfourd.*;
            for f = 1:length(this.hasEp2dDefaultMcf)
                try
                    niiMask = NIfTI(this.maskFile2(f));
                    niiHo   = NIfTI(this.hoFile(f));
                    niiHoSz = niiHo.size;
                    niiHoSz = niiHoSz(1:3);
                    if (~all(niiMask.size == niiHoSz(1:3)))
                        inc = [inc this.hasEp2dDefaultMcf{f}];
                        fprintf('%s:  size(mask) = %s, size(ho) = %s\n', this.hasEp2dDefaultMcf{f}, num2str(niiMask.size), num2str(niiHo.size));
                    else
                        this.viewMaskTransformation(f);
                    end
                catch ME %#ok<NASGU>
                    fprintf('SessionsCatalog.checkConsistencies encountered problems in %s\n', this.hasEp2dDefaultMcf{f});
                end
            end
        end
        function        rebuildMask(this, f)
            cd(fullfile(this.npPath, this.hasEp2dDefaultMcf{f}, 'fsl'));
            flirt = '/usr/local/fsl/bin/flirt';
            args  = '-bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 6  -interp trilinear';
            args2 = '-paddingsize 0.0 -interp nearestneighbour';
            system(sprintf('%s -in bt1_default.nii.gz -ref ho_meanvol_default.nii.gz -out bt1_default_on_ho_meanvol_default.nii.gz -omat bt1_default_on_ho_meanvol_default.mat %s', flirt, args));
            system(sprintf('%s -in bt1_default_mask.nii.gz  -applyxfm -init bt1_default_on_ho_meanvol_default.mat -out bt1_default_mask_on_ho_meanvol_default.nii.gz %s -ref ho_meanvol_default.nii.gz', flirt, args2));
            this.viewMaskTransformation(f);
        end
        function        viewMaskTransformation(this, f)
            cd(fullfile(this.npPath, this.hasEp2dDefaultMcf{f}, 'fsl'));
            system('fslview bt1_default_on_ho_meanvol_default ho_meanvol_default &');
        end
        function m    = missingHoMasked(this)
            m = this.findSomethingMissing(@this.hoMasked);
        end
  		function this = SessionsCatalog(varargin) 
 			%% SESSIONSCATALOG 
 			%  Usage:  this = SessionsCatalog(npPath) 

            p = inputParser;
            addOptional(p, 'npPath', pwd, @(x) lexist(x, 'dir'));
            parse(p, varargin{:});
            
            this.npPath = p.Results.npPath;
            cd(this.npPath);
            this = this.init_hasEp2dDefaultMcf;
 		end 
    end 
    
    %% PRIVATE
    
    methods (Access = 'private')
        function this = init_hasEp2dDefaultMcf(this)
            cd(this.npPath);
            dt = mlfourd.DirTool('mm0*');
            this.hasEp2dDefaultMcf = {};
            for t = 1:length(dt)
                if (lexist(fullfile(dt.dns{t}, this.FSL_FOLDER, 'ep2d_default_mcf.nii.gz'), 'file'))
                    this.hasEp2dDefaultMcf = [this.hasEp2dDefaultMcf dt.dns{t}];
                end
            end
        end
        function m    = findSomethingMissing(this, hfqfilename)
            m = {};
            for f = 1:length(this.hasEp2dDefaultMcf)
                if (~lexist(hfqfilename(f)))
                    m = [m this.hasEp2dDefaultMcf{f}];
                    fprintf('%s is missing %s\n', this.hasEp2dDefaultMcf{f}, hfqfilename(f));
                end
            end
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

