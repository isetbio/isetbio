classdef osLinear < outerSegment
% Linear subclass of the outersegment object
%
%   os = osLinear();
%
% Implements isomerizations (R*) to photocurrent (pA) using only cone
% linear temporal filters. The default values are those determined by
% the Angueyra and Rieke (2013, Nature Neuroscience).
%
% The current is calculated by convolving separate temporal filters
% for the L, M and S cones with the isomerization time course.
%
% JRG/HJ/BW, ISETBIO Team, 2016
    
    % Public, read-only properties.
    properties (GetAccess = public, SetAccess = private)
        % linear filters for LMS cones in columns
        lmsConeFilter;
    end
    
    properties (Access = private)
        % The properties in this part is used to store the state of the
        % outersegment so that we can do incremental computation
        
        absHistory = [];  % recent absorption history 
    end
    
    % Defined in separate files within the @osLinear directory
    methods (Access=public)
        [lmsFilters, meanCurrent] = generateBioPhysFilters(os, meanRate, varargin);
    end
    
    % Public methods
    methods
        % Constructor
        function obj = osLinear(varargin)
            % Initialize the parent class
            obj = obj@outerSegment();
        end
        
        % set function, see osLinearSet for details
        function obj = set(obj, varargin)
            osSet(obj, varargin{:});
        end
        
        % get function, see osLinearGet for details
        function val = get(obj, varargin)
            val = osGet(obj, varargin{:});
        end
        
        matchSensor(obj, varargin);
        
        function obj = compute(obj, pRate, coneType, varargin)
            % see osCompute for details
            obj = osCompute(obj, pRate, coneType, varargin{:});
        end
        
        function plot(obj, cMosaic)
            % see osPlot for details
            osPlot(obj, cMosaic);
        end
        
    end
    
    methods (Access = private)
        function newIRFs = generateLinearFilters(obj, meanRate, varargin)
            % To deprecate: Generates the linear temporal filters for the L, M and S
            % cones based on physiological data from Angueyra and Rieke
            % (Nat. Neuro., 2013).
            %
            %   newIRFs = generateLinearFilters([coneMosaic])
            %
            % Inputs: 
            %   obj      - osLinear class object
            %   meanRate - mean photon isomerization rate
            %
            % Outputs:
            %  newIRFs  - nx3 matrix with each row the temporal impulse
            %             response for the L, M and S cones, respectively.
            %
            % See also:
            %   osFitlerConesLinear
            %
            % JRG/HJ/BW, ISETBIO TEAM, 2016
            
            % Handle input
            p = inputParser;
            p.addRequired('obj', @(x) isa(x, 'osLinear'));
            p.addParameter('osType',0,@islogical);
            p.addRequired('meanRate', @isscalar);
            
            p.parse(obj, meanRate, varargin{:});          
                        
            osType = p.Results.osType; % peripheral (0) or foveal (1)
            
            dt = obj.timeStep;
            timeAxis = 0 : dt : 0.3;
            
            % parameters of cone temporal filters measured in Angueyra and
            % Rieke (Nat. Neuro., 2013).
            scFact = 1;      % scale factor
            tauR = 0.0216;   % rising phase time constant
            tauD = 0.0299;   % damping time constant
            tauP = 0.5311;   % period
            phi = 34.1814;   % phase
            
            % generate filter
            filter = scFact * (timeAxis/tauR).^3 ...
                ./ (1 + (timeAxis/tauR).^3) .* exp(-(timeAxis./tauD)) ...
                .* cos(2*pi*timeAxis/tauP+ deg2rad(phi));
            
            % Adjust gain of the cone impulse response for L, M, S cones
            %
            % Weber-Fechner relation:
            %  gain/gain_dark = ...
            %    1/(1+(bgIntensity/halfDesensitizationIntensitity))
            %            
            Io = 2250; % half-desensitizing background in R*/cone/sec
            Ib = [7131 6017 1973]; % R* per sec due to background adapting
            
            % Set gain parameters according to Angueyra & Rieke:
            Ib = Ib * meanRate / max(Ib);
            gain_dark = 0.22;              % from Juan's paper
            
            % scale IRF to reflect amplitude at chosen background
            % using Weber adaptation equation
            gain = gain_dark ./ (1 + Ib./Io);
            
            % The L, M and S IRFs are created by scaling Filter with the
            % scale factors
            newIRFs = filter(:) * gain * dt ./ max(filter);
            obj.lmsConeFilter = newIRFs;
        end
        
        
    end
end

