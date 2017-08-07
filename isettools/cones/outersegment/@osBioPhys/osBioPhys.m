classdef osBioPhys < outerSegment 
%osBioPhys  Create a biophysically based outersegment (os) object.
%
% Syntax:
%     os = osBioPhys;
%     os = osBioPhys('eccentricity',0);
%
% Description:
%     This class provides methods and parameters for converting a movie of
%     isomerizations (R*) to outer segment current (pA).
% 
%     Rieke and colleagues defined a set of difference equations as a
%     simulation of the phototransduction cascade. This object defines the
%     parameters and methods to transform the computed isomerizations (R*) in
%     the coneMosaic current.
%
%     If the noiseFlag property of the osLinear object is true, this method
%     adds noise to the current output signal.
%
%     The osBioPhys model is also the basis of how we find the linear filters
%     in the osLinear model, another subclass of outerSegment.
%
%     At present, we just have two sets of parameters, foveal and peripheral.  One day
%     we might try to handle eccentricity more finely.
%
% Optional key/value pairs:
%     'eccentricity'                Eccentricity in degrees. Determines parameters used.  Currently
%                                   we just have foveal and peripheral parameters, and somewhat
%                                   arbitrarily set the cuttoff at 2 degrees.  Default is 10 degrees.
%                                   
% References:
%     http://isetbio.org/cones/adaptation%20model%20-%20rieke.pdf
%     https://github.com/isetbio/isetbio/wiki/Cone-Adaptation


% JRG/HJ/BW, ISETBIO Team, 2016
%
% 08/05/17  dhb   Add eccentricityDegs key/value pair, and comments about it.
%                 This had been suggested previously in the comments as required.
%                 I did it in a backwards compatiable fashion.
%           dhb   Deprecate 'osType' and change 'eccentricityDegs' to 'eccentricity'.

    properties(Access = private)
        %state  Biophysics parameter state
        state; 
        
        %fovealPeripheralCutoffDegs  Eccentricity in degrees beyond which we switch from foveal to peripheral parameters
        % The 10 degrees value comes from Fred Rieke. Soon we will have
        % intermediate eccentricity dynamics (3-5 degrees) so this will
        % change.
        fovealPeripheralCutoffDegs = 10;
    end
    
    properties(SetAccess = protected, GetAccess = public)
        %model  Structure with biophysical model parameters
        model;        
    end  
   
    methods
        % Constructor
        function obj = osBioPhys(varargin)
            % Initialize the osBioPhys
           
            % Parse input
            p = inputParser;
            addParameter(p,'eccentricity',10,@isnumeric);
            p.parse(varargin{:});
            eccentricityDegs = p.Results.eccentricity; 
            
            % If eccentricity is greater than foveal/peripheral cutoff, use peripheral parameters
            if (eccentricityDegs > obj.fovealPeripheralCutoffDegs) 
                    % Peripheral parameters
                    obj.model.sigma = 22;       % rhodopsin activity decay rate (1/sec) - default 22
                    obj.model.phi = 22;         % phosphodiesterase activity decay rate (1/sec) - default 22
                    obj.model.eta = 2000;	    % phosphodiesterase activation rate constant (1/sec) - default 2000
                    obj.model.gdark = 20.5;     % concentration of cGMP in darkness - default 20.5
                    obj.model.k = 0.02;         % constant relating cGMP to current - default 0.02
                    obj.model.h = 3;            % cooperativity for cGMP->current - default 3
                    obj.model.cdark = 1;        % dark calcium concentration - default 1
                    obj.model.beta = 9;	        % rate constant for calcium removal in 1/sec - default 9
                    obj.model.betaSlow = 0.4;   % rate constant for slow calcium modulation of channels - default 0.4
                    obj.model.n = 4;  	        % cooperativity for cyclase, hill coef - default 4
                    obj.model.kGc = 0.5;        % hill affinity for cyclase - default 0.5
                    obj.model.OpsinGain = 10;   % So stimulus can be in R*/sec (this is rate of increase in opsin activity per R*/sec) - default 10
                    
            % Othewise use foveal parameters
            else
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
                    obj.model.OpsinGain = 12;   % So stimulus can be in R*/sec (this is rate of increase in opsin activity per R*/sec) - default 10
            end
                     
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
        
        function uData = plot(obj, sensor, varargin)
            % see osPlot for details
            uData = osPlot(obj, sensor, varargin{:});
        end
        
        function val = timeAxis(obj)
            % The temporal samples for the lms filters            
            val = ((1:size(linearFilters(obj.os,obj),1)) - 1) * obj.timeStep;
        end
        
    end
end
