classdef (Sealed) Fung2013Kit < handle & mlkinetics.InputFunctionKit
    %% line1
    %  line2
    %  
    %  Created 09-Jun-2022 14:20:06 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlaif/src/+mlaif.
    %  Developed on Matlab 9.12.0.1956245 (R2022a) Update 2 for MACI64.  Copyright 2022 John J. Lee.
    
    methods (Static)
        function this = instance(varargin)
            %% INSTANCE
            %  @param optional qualifier is char \in {'initialize' ''}
            
            ip = inputParser;
            addOptional(ip, 'qualifier', '', @ischar)
            parse(ip, varargin{:})
            
            persistent uniqueInstance
            if (strcmp(ip.Results.qualifier, 'initialize'))
                uniqueInstance = [];
            end          
            if (isempty(uniqueInstance))
                this = mlaif.Fung2013Kit();
                uniqueInstance = this;
            else
                this = uniqueInstance;
            end
        end
    end 

    %% PRIVATE

    methods (Access = private)
        function this = Fung2013Kit()
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
