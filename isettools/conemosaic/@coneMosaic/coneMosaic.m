classdef coneMosaic < handle
    % Create a cone mosaic class
    %
    %   coneM =  coneMosaic('cone',cone,'os','os);
    %
    % The cone mosaic defines an array of cones.  The individual cones have
    % absorption properties defined by coneM.cone.  The computation from
    % absorptions to cone voltages are defined by coneM.os
    %
    % JRG/BW ISETBIO Team, 2016
    
    
    properties (SetAccess = private, GetAccess = public)
        % These are the linear filters generated below via filterKernel.
        name = 'human-0';
        absorptions;      % The spatial array of cone absorptions
        
        cone;             % Properties of the individual cones
        os;
        
        wave;
        pattern;
        color;
        cfa;
        integrationTime;
        noiseFlag;
        
    end
    
    % Private properties. Only methods of the parent class can set these
    properties(Access = private)
    end
    
    % Public methods
    methods
        
        % Constructor
        function obj = coneMosaic(varargin)
            % Initialize the cone mosaic class
            p = inputParser;
            
            p.addParameter('cone',coneCreate,@isstruct);
            vFunc = @(x) (isa(x,'outerSegment'));
            p.addParameter('os',osCreate('linear'),vFunc);
            
            p.parse(varargin{:});
            obj.cone = p.Results.cone;
            obj.os   = p.Results.os;
            
            % Initialize the mosaic properties
            % obj.matchSensor(varargin{:});
            
        end
        
        % set function, see osLinearSet for details
        function obj = set(obj, varargin)
        end
        
        % get function, see osLinearGet for details
        function val = get(obj, varargin)
        end
        
    end
    
    % Methods that must only be implemented (Abstract in parent class).
    methods (Access=public)
                
        function obj = compute(obj, oi, varargin)
            % Compute the pattern of cone absorptions and typically the
            % photocurrent.  If you don't want the current, say so.
            p = inputParser;
            p.addParameter('current',true,@islogical)
            p.parse(varargin{:});
            current = p.Results.current;
            
            obj.computeAbsorptions(oi);
            
            % If we don't want the current, we can say so
            if current, obj.computeCurrent; end
        end
        
        function plot(obj, sensor)
        end
    end
    
    % Methods may be called by the subclasses, but are otherwise private
    methods (Access = protected)
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
    end
    
end

%
%
% switch sensorName
%     case {'default', 'human'}
%         % s = sensorCreate('human', [coneP], [retinalPos], [whichEye]);
%         % retinalPos should be 1x2 vector containing eccentricity (deg) and
%         % polar angle (deg)
%         if notDefined('coneP'), coneP = coneCreate; end
%         if ~isempty(varargin)
%             retPos = varargin{1};
%             if isscalar(retPos), retPos = [retPos, 0]; end
%         else
%             retPos = [0, 0];
%         end
%         if length(varargin)>1, whichEye = varargin{2};
%         else whichEye = []; end
%
%         eccMM = 2*tand(retPos(1)/2) * 17; % assuming focal length of 17 mm
%
%         % Assign key fields
%         wave = coneGet(coneP,'wave');
%         hPixel = pixelCreate('human', wave);
%
%         % Adjust pixel gap by retinal position
%         coneD = coneDensity(eccMM, retPos(2), whichEye);
%         coneSz = sqrt(1/coneD) * 1e-3; % avg cone size with gap in meters
%
%         % Adjust pixel gap size according to retinal position
%         wGap = coneSz - pixelGet(hPixel, 'width');
%         hGap = coneSz - pixelGet(hPixel, 'height');
%         assert(wGap>=0 && hGap>=0, 'gap should be non-negative');
%
%         hPixel = pixelSet(hPixel, 'width gap', wGap);
%         hPixel = pixelSet(hPixel, 'height gap', hGap);
%
%         % Add the default human pixel to the sensor
%         sensor.spectrum.wave = wave;
%         sensor = sensorSet(sensor, 'pixel', hPixel);
%         sensor = sensorSet(sensor, 'size', [72 88]);
%
%         % Add the default lens structure
%         lens = lensCreate([], wave);
%         sensor = sensorSet(sensor, 'human lens', lens);
%
%         % Add the default macular structure
%         macular = macularCreate(macularDensity(retPos(1)), wave);
%         sensor = sensorSet(sensor, 'human macular', macular);
%
%         % Build up a human cone mosaic.
%         sensor = sensorCreateConeMosaic(sensor, coneP);
%



