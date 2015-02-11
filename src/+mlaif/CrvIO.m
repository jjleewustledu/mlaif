classdef CrvIO < mlio.AbstractIO
	%% CRVIO is a concrete class for filesystem I/O of ASCII/Unicode text
    %    
	%  $Revision$
 	%  was created $Date$
 	%  by $Author$, 
 	%  last modified $LastChangedDate$
 	%  and checked into repository $URL$, 
 	%  developed on Matlab 8.1.0.604 (R2013a)
 	%  $Id$

    properties (Constant)
        FILETYPE     = 'ECAT EXACT HR+ arterial line';
        FILETYPE_EXT = '.crv';
    end
    
    properties (Dependent)
        array
        contents
        counts
        descrip
        dt
        header
        time
    end
    
	methods (Static)
        function this  = load(fn) 
            assert(lexist(fn, 'file'));
            [pth, fp, fext] = fileparts(fn); 
            if (strcmp(mlaif.CrvIO.FILETYPE_EXT, fext) || ...
                isempty(fext))
                this = mlaif.CrvIO.loadText(fn); 
                this.fqfilename = fullfile(pth, [fp fext]);
                return 
            end
            error('mlaif:unsupportedParam', 'CrvIO.load does not support file-extension .%s', fext);
        end
        function ca    = textfileToCell(fqfn, eol)  %#ok<INUSD>
            if (~exist('eol','var'))
                fget = @fgetl;
            else
                fget = @fgets;
            end
            ca = {[]};
            try
                fid = fopen(fqfn);
                i   = 1;
                while 1
                    tline = fget(fid);
                    if ~ischar(tline), break, end
                    ca{i} = tline;
                    i     = i + 1;
                end
                fclose(fid);
                assert(~isempty(ca) && ~isempty(ca{1}))
            catch ME
                fprintf('mlaif.CrvIO.textfileToCell:  exception thrown while reading \n\t%s\n\tME.identifier->%s', fqfn, ME.identifier);
            end
        end
    end
    
    methods %% GET
        function a     = get.array(this)
            a = this.contents(3:end);
            a = cellfun(@str2num, a, 'UniformOutput', false)';
            a = cell2mat(a);
        end
        function c     = get.contents(this)
            assert(~isempty(this.contents_));
            c = this.contents_;
        end
        function c     = get.counts(this)
            c = this.array(:,2);
        end
        function d     = get.descrip(this)
            d = sprintf('%s read %s on %s', class(this), this.fqfilename, datestr(now));
        end
        function t     = get.dt(this)            
            rgxnames = regexp(this.header, this.binWidthRegexp_, 'names');
            t = str2double(rgxnames.binwidth);
        end
        function h     = get.header(this)
            h = [this.contents{1} ' ' this.contents{2}];
        end
        function t     = get.time(this)
            t = this.dt * this.array(:,1);
        end
    end
    
    methods
        function ch    = char(this)
            ch = strjoin(this.contents, '\n');
        end
        function         plot(this)
            figure;
            plot(this.time, this.counts);
        end
        function         save(~)
        end
    end
    
    %% PRIVATE
    
    properties (Access = 'private')
        contents_
        binWidthRegexp_ = 'BinWidth=(?<binwidth>\d+.\d+)\s+seconds\s+\d+\s+\d+';
    end
    
    methods (Static, Access = 'private')
        function this = loadText(fn)
            import mlaif.*;
            this = CrvIO;
            this.contents_ = CrvIO.textfileToCell(fn);
        end
    end
    

	%  Created with Newcl by John J. Lee after newfcn by Frank Gonzalez-Morphy 
end

