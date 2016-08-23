classdef osBioPhys < outerSegment 
% BioPhys subclass of the outersegment object
% 
%       os = osBioPhys();
% 
% Converts isomerizations (R*) to outer segment current (pA). The
% difference equation model by Rieke implements a biophysical
% simulation of the phototransduction cascade. If the noiseFlag
% property of the osLinear object is set to 1, this method will add noise
% to the current output signal.
%
% Reference:
%   http://isetbio.org/cones/adaptation%20model%20-%20rieke.pdf
%   https://github.com/isetbio/isetbio/wiki/Cone-Adaptation
% 
% JRG/HJ/BW, ISETBIO Team, 2016

    properties(Access = private)
        state;   % biophysics parameter state
    end
    
    properties(SetAccess = protected, GetAccess = public)
        model;
        
%         sigma;   % rhodopsin activity decay rate (1/sec) - default 22
%         phi;     % phosphodiesterase activity decay rate (1/sec) - default 22
%         eta;     % phosphodiesterase activation rate constant (1/sec) - default 2000
%         gdark;   % concentration of cGMP in darkness - default 20.5
%         k;       % constant relating cGMP to current - default 0.02
%         h;       % cooperativity for cGMP->current - default 3
%         cdark;   % dark calcium concentration - default 1
%         beta;    % rate constant for calcium removal in 1/sec - default 9
%         betaSlow;% rate constant for slow calcium modulation of channels - default 0.4
%         n;       % cooperativity for cyclase, hill coef - default 4
%         kGc;     % hill affinity for cyclase - default 0.5
%         OpsinGain; % so stimulus can be in R*/sec (this is rate of increase in opsin activity per R*/sec) - default 10
%         
%         q;      % constant accoutning for fraction of the current carried by calcium and the outer segment volume
%         smax;   % max cGMP creation rate
%         
%         bgCur;  % background current
%         
%         opsin   % photopigment activity
%         PDE     % cGMP hydrolysis rate by phoshpodiesterase 
%         Ca      % calcium concentration
%         Ca_slow % lowpass-filtered calcium concentration
%         st      % cGMP creation rate
%         cGMP    % cGMP concentration
    end
    

    
    methods
        
        function obj = osBioPhys(varargin)
            
            p = inputParser;
            addParameter(p,'osType',0,@islogical);
            p.parse(varargin{:});
            
            osType = p.Results.osType; % peripheral (0) or foveal (1)
            
            switch osType
                
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
        
        % @JRG need to make new object version of osAddNoise
        [adaptedCur, params] = osAddNoise(curNF, params);
    end
    
    methods (Access=public)        
        function obj = compute(obj, sensor, varargin)
            % see osCompute for details
            obj = osCompute(obj, sensor, varargin{:});
        end
        
        function plot(obj, sensor, varargin)
            % see osPlot for details
            osPlot(obj, sensor, varargin{:});
        end
    end
end
