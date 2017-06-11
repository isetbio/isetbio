classdef osBioPhys < outerSegment 
% The parameters and methods to convert isomerizations (R*) to outer
% segment current (pA).
% 
%       os = osBioPhys();
% 
% Rieke and colleagues defined a set of difference equations as a
% simulation of the phototransduction cascade. This object defines the
% parameters and methods to transform the computed isomerizations (R*) in
% the coneMosaic current.
%
% If the noiseFlag property of the osLinear object is true, this method
% adds noise to the current output signal.
%
% The osBioPhys model is also the basis of how we find the linear filters
% in the osLinear model, another subclass of outerSegment.
%
% Reference:
%   http://isetbio.org/cones/adaptation%20model%20-%20rieke.pdf
%   https://github.com/isetbio/isetbio/wiki/Cone-Adaptation
% 
% TODO:
%  See constructor comments about osType (BW)
%
% JRG/HJ/BW, ISETBIO Team, 2016

    properties(Access = private)
        state;   % biophysics parameter state
    end
    
    properties(SetAccess = protected, GetAccess = public)
        model;        
    end  
   
    methods
        
        function obj = osBioPhys(varargin)
            
            p = inputParser;
            addParameter(p,'osType',0,@islogical);
            p.parse(varargin{:});
            
            % We should rename the osType parameter and make eccentricity a
            % number. Everything less than X should be foveal and greater
            % than X should be peripheral.
            eccentricity = p.Results.osType; % peripheral (0) or foveal (1)
            
            switch eccentricity
                
                case 0 % peripheral
                    % Peripheral parameters
                    obj.model.sigma = 22;  % rhodopsin activity decay rate (1/sec) - default 22
                    obj.model.phi = 22;     % phosphodiesterase activity decay rate (1/sec) - default 22
                    obj.model.eta = 2000;	  % phosphodiesterase activation rate constant (1/sec) - default 2000
                    obj.model.gdark = 20.5; % concentration of cGMP in darkness - default 20.5
                    obj.model.k = 0.02;     % constant relating cGMP to current - default 0.02
                    obj.model.h = 3;       % cooperativity for cGMP->current - default 3
                    obj.model.cdark = 1;  % dark calcium concentration - default 1
                    obj.model.beta = 9;	  % rate constant for calcium removal in 1/sec - default 9
                    obj.model.betaSlow = 0.4; % rate constant for slow calcium modulation of channels - default 0.4
                    obj.model.n = 4;  	  % cooperativity for cyclase, hill coef - default 4
                    obj.model.kGc = 0.5;   % hill affinity for cyclase - default 0.5
                    obj.model.OpsinGain = 10; % so stimulus can be in R*/sec (this is rate of increase in opsin activity per R*/sec) - default 10
          
                case 1 % foveal
                    % Foveal parameters
                    obj.model.sigma = 10;       % rhodopsin activity decay rate (1/sec) - default 22
                    obj.model.phi   = 22;       % phosphodiesterase activity decay rate (1/sec) - default 22
                    obj.model.eta   = 700;      % phosphodiesterase activation rate constant (1/sec) - default 2000
                    obj.model.gdark = 20.5;     % concentration of cGMP in darkness - default 20.5
                    obj.model.k     = 0.02;     % constant relating cGMP to current - default 0.02
                    obj.model.h     = 3;        % cooperativity for cGMP->current - default 3
                    obj.model.cdark = 1;        % dark calcium concentration - default 1
                    obj.model.beta  = 5;        % rate constant for calcium removal in 1/sec - default 9
                    obj.model.betaSlow = 0.4;   % rate constant for slow calcium modulation of channels - default 0.4
                    obj.model.n     = 4;        % cooperativity for cyclase, hill coef - default 4
                    obj.model.kGc   = 0.5;      % hill affinity for cyclase - default 0.5
                    obj.model.OpsinGain = 12;   % so stimulus can be in R*/sec (this is rate of increase in opsin activity per R*/sec) - default 10
            end
            
            % % Derived properties
            % obj.model.q    = 2 * obj.model.beta * obj.model.cdark / (obj.model.k * obj.model.gdark^obj.model.h);
            % obj.model.smax = obj.model.eta/obj.model.phi * obj.model.gdark * (1 + (obj.model.cdark / obj.model.kGc)^obj.model.n);
                     
        end
        
        function obj = set(obj, varargin)
            % set function, see osBioPhysSet for details
            osSet(obj, varargin{:});
        end
        
        
        function val = get(obj, varargin)
            % get function, see osBioPhysGet for details
            val = osGet(obj, varargin{:});
        end
        
        state = osAdaptSteadyState(obj, bgR, varargin);
        
        [adaptedData, state] = osAdaptTemporal(pRate,obj);
        
    end
    
    methods (Access=public)        
        function obj = compute(obj, varargin)
            % see osCompute for details
            obj = osCompute(obj, varargin{:});
        end
        
        function plot(obj, sensor, varargin)
            % see osPlot for details
            osPlot(obj, sensor, varargin{:});
        end
        
        function val = timeAxis(obj)
            % The temporal samples for the lms filters
%             if isempty(obj.lmsConeFilter)
%                 warning('No lms impulse response functions computed');
%             else
                % Time axis is the length of the filters multiplied by the
                % time step
                val = ((1:size(linearFilters(obj.os,obj),1)) - 1) * obj.timeStep;
%             end
        end
        
    end
end
