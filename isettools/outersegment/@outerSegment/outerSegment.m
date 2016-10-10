classdef outerSegment < handle
% The outerSegment parent class for modeling the responses of the cone
% outer segment.
% 
% The subclasses osLinear and osBioPhys calculate outer segment current
% responses using either a linear temporal filter model or a biophysical
% difference equations model (both based on physiology from Fred Rieke's
% lab, and both models determined by Fred Rieke).
%
% osLinear implements isomerizations (R*) to photocurrent (pA) using only
% cone linear temporal filters. The default values are those determined by
% the Angueyra and Rieke (2013, Nature Neuroscience).
% 
% osBioPhys converts isomerizations (R*) to outer segment current (pA). The
% difference equation model by Rieke is applied here. If the noiseFlag
% property of the osLinear object is set to 1, this method will add noise
% to the current output signal.
% 
% Reference for osBioPhys model:
%   http://isetbio.org/cones/adaptation%20model%20-%20rieke.pdf
%   https://github.com/isetbio/isetbio/wiki/Cone-Adaptation
%
% The class also allows the addition of noise based on physiological
% measurements of cone responses to a neutral gray background (physiology
% and model also by Fred Rieke).
%
% See also subclasses:
%       osLinear.m, osBioPhys.m
% 
% JRG, NC, DHB, 8/2015
% 
    % Public read/write properties
    properties
    end
    
    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)  

    end
    
    % Protected properties; Methods of the parent class and all of its
    % subclasses can set these.
    %
    properties
        noiseFlag;            % determines whether noise is added
        patchSize;            % spacing between cones (width) in um
    end
    
    properties (SetAccess = protected)
        coneCurrentSignal;    % output signal in pA
    end
    
    properties (SetObservable, AbortSet)
        timeStep;             % sampling interval in sec
    end
    
    % Public methods
    methods       
        function obj = outerSegment(varargin)
            obj.noiseFlag = 0;            
            obj.coneCurrentSignal = [];
            obj.timeStep = 1e-3;
        end
        
        % see osSet in @osLinear and @osBioPhys for details
        function obj = set(obj, varargin)
            osSet(obj, varargin{:});
        end
        
        % see osGet in @osLinear and @osBioPhys for details
        function val = get(obj, varargin)
           val = osGet(obj, varargin{:});
        end
        
    end
    
    % Methods that must only be implemented in the subclasses.
    % If a subclass does not implement each and every of these methods
    % it cannot instantiate objects.
    methods (Abstract, Access=public)
        % see osLinearCompute, osBioPhysCompute
        compute(obj, pRate, coneType, varargin);
        % see osLinearPlot, osBioPhysPlot
        plot(obj, plotType);
    end
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)

    end
    
    methods (Static)
        resampledPhotocurrents = resample(photocurrents, originalTimeAxis, resampledTimeAxis);
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        initialize(obj);
    end
    
end
