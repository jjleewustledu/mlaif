classdef LaifPET <  mlaif.AbstractAifProblem 
	%% LAIFPET  
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
        HOME = '/Volumes/PassportStudio2/cvl/np755'
        PLOT_4 = true    
    end
       
    properties
        showPlots = true;
        studyId        
        xLabel    = 'time/s'
        yLabel    = 'counts'
        
        aifPET
        parmax
        avpar
        estimatedAif
        estimatedWb
    end
    
    properties (Dependent)
        baseTitle
        aifConc
    end
    
    methods %% GET
        function bt = get.baseTitle(this)
            bt = sprintf('PET AIF, whole-brain from %s', this.studyId);
        end
        function c  = get.aifConc(this)   
            assert(~isempty(this.aifPET));
            c = this.aifPET.tracerConcentrations;
        end
    end
    
    methods (Static)
        function bps = runSessions(varargin)
            cd(mlaif.LaifPET.HOME);
            p = inputParser;
            addOptional(p, 'tInitial', 1, @isnumeric);
            parse(p, varargin{:});
            tinit = p.Results.tInitial;
            
            load('bayesianPets_2014nov11.mat');
            bps = cell(1,length(sessStrings)); %#ok<USENS>
            for t = tinit:length(sessStrings)
                try
                    bps{t} = mlaif.LaifPET.runSession(sessStrings{t}, fracaif(t), t01(t), toffset(t)); 
                catch ME %#ok<NASGU>
                    warning('mlaif:sessionException', 'LaifPET.runSessions encountered errors in %s', sessStrings{t});
                end
            end
            diary off
        end
        function LaifPET  = runSession(sessStr, fracaif, t01, toffset)
            cd(mlaif.LaifPET.HOME);
            bayesDir = fullfile(sessStr, 'bayesian_pet', '');
            crvDir   = fullfile(sessStr, 'ECAT_EXACT', 'pet', '');
            assert(lexist(bayesDir, 'dir'), 'LaifPET.runSession could not find dir %s', bayesDir);
            assert(lexist(crvDir, 'dir'),   'LaifPET.runSession could not find dir %s', crvDir);
            
            pwd0 = pwd;            
            try
                import mlfourd.*;
                cd(bayesDir);                
                diary('bayesianPet_2014no11.log');
                dt = DirTool(fullfile(crvDir, 'p*ho*.crv'));
                LaifPET = mlaif.LaifPET(sessStr, strtok(dt.fns{1}, '.'));
                LaifPET = LaifPET.estimateParameters(fracaif, t01, toffset);
                filename = fullfile(bayesDir, ['bayesianPet_' str2mmnum(bayesDir) '_' str2pnum(bayesDir) '.mat']);
                save(filename, 'LaifPET');
                fprintf('LaifPET.runSession saved %s\n', filename);
                diary off
                cd(pwd0);
            catch ME
                handexcept(ME);
            end
            
        end
    end
    
	methods 		
        function this = estimateParameters(this, fracaif, t01, toffset)
            %% ESTIMATEPARAMETERS manages Bayes PETMR processing
            %  Usage:  data_object = this.estimateParameters
            %          ^ struct

            pmap = containers.Map;
            ncnt0 = max(this.aifPET.tracerConcentrations);
            pmap('c1') = struct('fixed', 0, 'min', -0.2, 'mean', -0.0125, 'max', 0);
            pmap('c2') = struct('fixed', 0, 'min',  0,   'mean', 15.80, 'max', 20);
            pmap('c3') = struct('fixed', 0, 'min',  0, 'mean', 1.044, 'max', 3);
            pmap('c4') = struct('fixed', 0, 'min',  0, 'mean', 0.9224, 'max', 5);
            pmap('fracSS') = struct('fixed', 0, 'min', 0.01, 'mean', 0.0120, 'max', 0.05);
            pmap('fracaif') = struct('fixed', 1, 'min', 0.0015, 'mean', fracaif, 'max', 0.005);
            pmap('k') = struct('fixed', 0, 'min', 0, 'mean', 0.1322, 'max', 0.14); % coefficient for PET CBF
            pmap('kappa') = struct('fixed', 0, 'min', 1, 'mean', 7.560, 'max', 11); % MR->PET time constant
            pmap('ncnt') = struct('fixed',0 , 'min', 0.4e4, 'mean', ncnt0, 'max', 2e4);            
            pmap('t01') = struct('fixed', 1, 'min', 15, 'mean', t01, 'max', 45);
            pmap('t02') = struct('fixed', 0, 'min', 0, 'mean', 47 , 'max', 100); % a duration
            pmap('toffset') = struct('fixed', 1, 'min', 10, 'mean', toffset, 'max', 26);         
            
            %% Call the MCMC for processing
            
            import mlbayesian.*;
            this.paramsManager = BayesianParameters(pmap);            
            this.mcmc          = MCMC(this, this.dependentData, this.paramsManager);
            [this.parmax,this.avpar,this.mcmc] = this.mcmc.runMcmc;   
            
            %% get final plots            
            
            [this.estimatedAif, this.estimatedWb] = this.estimateData;   
            this.plotEstimate;
            this.plotFour;
        end
        function        plotOri(this)
            if (~this.PLOT_ORI)
                return; end
            figure
            plot(this.independentData, this.aifConc / max(this.aifConc), ...
                 this.independentData, this.dependentData  / max(this.dependentData),  'b--')
            title(['PET AIF, whole-brain from ', this.studyId])
            xlabel('time/s')
            ylabel(sprintf('tracer concentration (%g, %g)', ...
                   max(this.aifConc), max(this.dependentData)))
        end
        function        plotEstimate(this)
            if (~this.PLOT_ESTIMATE)
                return; end
            figure
            plot(this.independentData, this.aifPET.tracerConcentrations, ...
                 this.independentData, this.estimatedAif, 'b--');
            title(sprintf('PET AIF from %s, Bayesian estimated AIF', this.studyId))
            xlabel('time/s')
            ylabel(sprintf('tracer concentration (%g, %g)', ...
                   max(this.aifPET.tracerConcentrations), max(this.estimatedAif)))
        end
        function        plotFour(this)
            if (~this.PLOT_4)
                return; end
            figure
            plot(this.independentData, this.aifConc / max(this.aifConc), ...
                 this.independentData, this.dependentData  / max(this.dependentData),                    'b--', ...
                 this.independentData, this.estimatedAif / max(this.estimatedAif), 'g-.', ...
                 this.independentData, this.estimatedWb  /  max(this.dependentData),  'r:');
            title(sprintf('PET AIF & whole-brain from %s; AIF & whole-brain estimated', this.studyId))
            xlabel('time/s')
            ylabel(sprintf('tracer concentration (%g, %g, %g, %g)', ...
                   max(this.aifConc), max(this.dependentData), max(this.estimatedAif), max(this.estimatedWb)))   
        end
        function sse = sumSquaredErrors(this, p)
            p = num2cell(p);
            [~,wbtac] = this.aifPET.doublePassFast(p{:});
            sse = norm(this.dependentData  - wbtac)^2  / norm(this.dependentData)^2;
        end
        
        function [aifConc,wb] = estimateData(this)
            [aifConc,wb] = this.estimateDataFast( ....
                this.finalParams('c1'),  this.finalParams('c2'),    this.finalParams('c3'),     this.finalParams('c4'), ...
                this.finalParams('fracSS'), this.finalParams('fracaif'), ...
                this.finalParams('k'),   this.finalParams('kappa'), this.finalParams('ncnt'), ...
                this.finalParams('t01'), this.finalParams('t02'),   this.finalParams('toffset'));
        end
        function [aifConc,wb] = estimateDataFast(this, ...
                c1, c2, c3, c4, fracSS, fracaif, k, kappa, ncnt, t01, t02, toffset)
            [aifConc,wb] = this.aifPET.doublePassFast(c1, c2, c3, c4, fracSS, fracaif, k, kappa, ncnt, t01, t02, toffset);
        end        
        
  		function this  = LaifPET(varargin)
 			%% LAIFPET 
 			%  Usage:  this = LaifPET()
 			
            p = inputParser;
            addRequired(p, 'sessionStr', @(x) assert(lexist(x, 'dir')));
            addRequired(p, 'studyId',    @(x) assert(lstrfind(x, 'p')));
            parse(p, varargin{:});
            sessStr      = p.Results.sessionStr;
            this.studyId = p.Results.studyId;
            
            perfDir  = fullfile(sessStr, 'perfusion_4dfp', '');
            hoDir    = fullfile(sessStr, 'ECAT_EXACT', '962_4dfp', '');
            crvDir   = fullfile(sessStr, 'ECAT_EXACT', 'pet', '');
            assert(lexist(perfDir, 'dir'));
            assert(lexist(hoDir, 'dir'));
            assert(lexist(crvDir, 'dir'));
            
            import mlaif.*;           
            this.aifPET          = AifPET.load( ...
                                   fullfile(crvDir, [this.studyId '.crv']), ...
                                   fullfile(perfDir, 'perfusion_4dfp.log'));
            this.independentData = this.aifPET.timeInterpolants;
            wbPET                = WholeBrainPET.load( ...
                                   fullfile(hoDir,  [this.studyId '_masked']));
            this.dependentData   = wbPET.tracerConcentrations;
 		end 
    end

	%  Created with NewClassStrategy by John J. Lee, after newfcn by Frank Gonzalez-Morphy 
end

