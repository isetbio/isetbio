classdef osLinear < outerSegment 
% @osLinear: Subclass of @outerSegment object
% 
% This subclass implements the cone linear filters determined by the
% experiments found in Angueyra and Rieke (2013) to convert cone 
% isomerizations (R*) to current (pA).
% 
% % Example code
% nSamples = 2000;        % 2000 samples
% timeStep = 1e-4;        % time step
% flashIntens = 50000;    % flash intensity in R*/cone/sec (maintained for 1 bin only)
% sensor = sensorCreate('human');
% sensor = sensorSet(sensor, 'size', [1 1]); % only 1 cone
% sensor = sensorSet(sensor, 'time interval', timeStep); 
% stimulus = zeros(nSamples, 1); stimulus(1) = flashIntens;
% stimulus = reshape(stimulus, [1 1 nSamples]);
% sensor = sensorSet(sensor, 'photon rate', stimulus);
% noiseFlag = 0; adaptedOS = osLinear('noiseFlag', noiseFlag);
% sensor = sensorSet(sensor,'adaptation offset',params.bgVolts);
% sensor = sensorSet(sensor,'cone type', 2); % set s cone
% adaptedOS = osLinearCompute(adaptedOS, sensor);
% osAdaptedCur = osLinearGet(adaptedOS, 'ConeCurrentSignal');
% osAdaptedCur = osAdaptedCur - osAdaptedCur(:, :, nSamples);
% osLinearPlot(adaptedOS, sensor);

    % Public, read-only properties.
    properties (SetAccess = private, GetAccess = public)
        % These are the linear filters generated below via filterKernel.
        sConeFilter;
        mConeFilter;
        lConeFilter;
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
    
    % Public methods
    methods
        
        % Constructor
        function obj = osLinear(varargin)
            % Initialize the parent class
            obj = obj@outerSegment();
            
            % Initialize ourselves
            obj.initialize();
            
            % parse the varargin
            for k = 1:2:numel(varargin)
                obj.(varargin{k}) = varargin{k+1};
            end
        end
        
        % set function, see osLinearSet for details
        function obj = set(obj, param, val, varargin)
            osLinearSet(obj, param, val, varargin);
        end
        
        % get function, see osLinearGet for details
        function val = get(obj, param, varargin)
           val = osLinearGet(obj, param, varargin);
        end
      
    end
    
    % Methods that must only be implemented (Abstract in parent class).
    methods (Access=public)
        function obj = compute(obj, sensor, param, varargin)
            % see osLinearCompute for details
            obj = osLinearCompute(obj, sensor, varargin); 
        end
        function plot(obj, sensor)
            % see osLinearPlot for details
            osLinearPlot(obj, sensor);
        end
    end    
    
    % Methods may be called by the subclasses, but are otherwise private 
    methods (Access = protected)
        % Generates the linear filters for L, M and S cones using data from
        % physiology by Angueyra and Rieke (2013)
        filterKernel(obj, varargin);
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        initialize(obj);
    end
    
end
