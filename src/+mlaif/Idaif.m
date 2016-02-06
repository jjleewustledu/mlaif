classdef Idaif 
	%% IDAIF  

	%  $Revision$
 	%  was created 17-Jan-2016 17:57:19
 	%  by jjlee,
 	%  last modified $LastChangedDate$
 	%  and checked into repository /Users/jjlee/Local/src/mlcvl/mlaif/src/+mlaif.
 	%% It was developed on Matlab 8.5.0.197613 (R2015a) for MACI64.
 	

	properties (Constant)
 		SUBJECTS = {'p*' 'NP995_*' 'HYGLY*'}
        VISITS = 'V*'
        SCANS = {'fdg*' 'ho*' 'oc*' 'oo*'}
        HDR = '*.4dfp.hdr'
 	end

	methods (Static)
        function [s,r] = subjects2niftigz()
            s = -1; r = '';
            subjects = mlsystem.DirTools(mlaif.Idaif.SUBJECTS);
            for u = 1:length(subjects)
                cd(subjects.dns{u});
                try                    
                    [s,r] = mlaif.Idaif.hdr2niftigz;
                    [s,r] = mlaif.Idaif.visits2nifitgz;
                catch ME
                    handwarning(ME);
                end
            end
        end
        function [s,r] = visits2niftigz()
            s = -1; r = '';
            visits = mlsystem.DirTool(mlaif.Idaif.VISITS);
            for v = 1:length(visits)
                cd(visits.dns{v});
                try
                    [s,r] = mlaif.Idaif.hdr2niftigz;
                    [s,r] = mlaif.Idaif.scans2nifitgz;
                catch ME
                    handwarning(ME);
                end
            end            
        end
        function [s,r] = scans2niftigz()
            s = -1; r = '';
            scanfolders = mlsystem.DirTools(mlaif.Idaif.SCANS);
            for d = 1:length(scanfolders.fqdns)
                try     
                    fprintf('scans2niftigz working in %s..........\n', scanfolders.fqdns{d});
                    cd(scanfolders.fqdns{d});
                    [s,r] = mlaif.Idaif.hdr2niftigz;
                    cd('pet_proc');
                    [s,r] = mlaif.Idaif.hdr2niftigz;
                catch ME
                    handwarning(ME);
                end
            end
        end
 		function [s,r] = hdr2niftigz()
            s = -1; r = '';
            hdrfiles = mlsystem.DirTool(mlaif.Idaif.HDR);
            for f = 1:length(hdrfiles.fqfns)
                try
                    fprintf('hdr2niftigz working on %s..........\n', hdrfiles.fqfns{f});
                    [s,r] = mlbash(sprintf('fslchfiletype NIFTI_GZ %s', hdrfiles.fqfns{f}));
                catch ME
                    handwarning(ME);
                end
            end
 		end
    end 
    
    methods 
        function this = Idaif()
        end
    end

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy
 end

