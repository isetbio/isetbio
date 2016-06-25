classdef coneMosaic < handle
    % Create a cone mosaic class
    %
    %   cMosaic =  coneMosaic('cone', cone, 'os', os);
    %
    % The cone mosaic defines an array of cones.  The individual cones have
    % absorption properties defined by cMosaic.cone.  The computation from
    % absorptions to photocurrent is defined by cMosaic.os
    %
    % HJ/JRG/BW ISETBIO Team, 2016
    
    properties  % public properties
        
        cone;             % Cone class object, contain single cone property
        macular;          % Macular class object
        os;               % Outersegment properties
        
        pattern;          % Pattern of K-LMS cones in the mosaick
        integrationTime;  % Cone temporal integration time in secs
        positions;        % Eye movement positions in number of cones.
                          % The length of this property controls number of
                          % frames to be computed
        noiseFlag;        % To control which noise is included
    end
    
    properties (GetAccess=public, SetAccess=private)
        absorptions;      % The spatial array of cone absorptions over time
        current;          % The spatial array of photocurrent over time
    end
    
    properties (Dependent)
        wave;       % Wavelength samples
        
        rows;       % number of rows in the cone mosaic
        cols;       % number of cols in the cone mosaic
        mosaicSize; % [rows, cols]
        
        width;      % width of cone mosaic in meters
        height;     % height of cone mosaic in meters
        
        qe;         % effective absorptance with macular pigment
    end
    
    % Public methods
    methods
        
        % Constructor
        function obj = coneMosaic(varargin)
            % Initialize the cone mosaic class
            %
            %   cMosaic =  coneMosaic('cone',cone,'os','os);
            %
            
            % parse input
            p = inputParser;
            
            p.addParameter('cone', Cones(),@(x) isa(x, 'Cones'));
            p.addParameter('macular', Macular(), @(x)isa(x, 'Macular'));
            p.addParameter('wave', 400:10:700, @isnumeric);
            p.addParameter('integrationTime', 0.05, @isscalar);
            p.addParameter('positions', [0 0], @isnumeric);
            p.addParameter('noiseFlag', 1, @isscalar);
            p.addParameter('pattern', [], @isnumeric);
            p.addParameter('size', [72 88], @isnumeric);
            p.addParameter('os', osCreate('linear'), ...
                @(x)(isa(x,'outerSegment')));
            
            p.parse(varargin{:});
            
            % set properties
            obj.cone = p.Results.cone;
            obj.macular = p.Results.macular;
            obj.os = p.Results.os;
            obj.wave = p.Results.wave;
            obj.integrationTime = p.Results.integrationTime;
            obj.noiseFlag = p.Results.noiseFlag;
            obj.positions = p.Results.positions;
            
            if isempty(p.Results.pattern)
                [~, obj.pattern] = humanConeMosaic(p.Results.size, ...
                    obj.cone.spatialDensity, obj.cone.width);
            else
                obj.pattern = p.Results.pattern;
            end
            
            % Initialize the mosaic properties
            % obj.matchSensor(varargin{:});
            
            % initialize listener
            % these listeners are used to make sure the wavelength samples
            % in obj.cone and obj.macular are always the same
            addlistener(obj.cone, 'wave', 'PostSet', @obj.setWave);
            addlistener(obj.macular, 'wave', 'PostSet', @obj.setWave);
        end
        
        
        % get methods for dependent variables
        function val = get.wave(obj)
            val = obj.cone.wave;
        end
        
        function val = get.rows(obj)  % number of rows
            val = size(obj.pattern, 1);
        end
        
        function val = get.cols(obj)  % number of cols
            val = size(obj.pattern, 2);
        end
        
        function val = get.width(obj)  % width of cone mosaic in meters
            val = obj.cols * obj.cone.width;
        end
        
        function val = get.height(obj)  % height of cone mosaic in meters
            val = obj.rows * obj.cone.height;
        end
        
        function val = get.qe(obj)
            % compute effective absorptance with macular pigments
            val = bsxfun(@times, obj.cone.absorptance, ...
                obj.macular.transmittance) * diag(obj.cone.peakEfficiency);
        end
        
        % set method for dependent variables
        function set.wave(obj, val)
            obj.cone.wave = val(:);
            obj.macular.wave = val(:);
        end
        
        function set.mosaicSize(obj, val)
            if any(val ~= obj.mosaicSize)
                [~, obj.pattern] = humanConeMosaic(val, ...
                    obj.cone.spatialDensity, obj.cone.width); 
            end
        end
        
        function set.rows(obj, val)
            obj.mosaicSize = [val obj.cols];
        end
        
        function set.cols(obj, val)
            obj.mosaicSize = [obj.rows val];
        end
        
        % set function, see osLinearSet for details
        % function obj = set(obj, varargin)
        % end
        
        % get function, see osLinearGet for details
        % function val = get(obj, varargin)
        %    val = obj;  % Place holder for what will happen!
        % end
        
    end
    
    % Methods that must only be implemented (Abstract in parent class).
    methods (Access=public)
        
        function obj = compute(obj, oi, varargin)
            % coneMosaic.plot()
            %
            % Compute the pattern of cone absorptions and typically the
            % photocurrent.  If you don't want the current, say so.
            p = inputParser;
            p.addParameter('current',true,@islogical)
            p.parse(varargin{:});
            current = p.Results.current;
            
            % Places the absorptions in the absorptions slot
            obj.computeAbsorptions(oi);
            
            % If we don't want the current, we can say so
            % Otherwise, place the current in the current slot
            if current, obj.computeCurrent; end
        end
        
        function plot(obj,varargin)
            % coneMosaic.plot()
            p = inputParser;
            % Do some plotting based on the input arguments.
        end
    end
    
    % Methods may be called by the subclasses, but are otherwise private
    methods (Access = protected)
    end
    
    % Methods that are totally private (subclasses cannot call these)
    methods (Access = private)
        % callback function for listeners
        function setWave(obj, src, ~)
            switch src.DefiningClass.Name
                case 'Cones'
                    obj.macular.wave = obj.cone.wave;
                case 'Macular'
                    obj.cone.wave = obj.macular.wave;
            end
        end
    end
    
end