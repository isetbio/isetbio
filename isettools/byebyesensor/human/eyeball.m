classdef eyeball < handle
    %% Class eyeball
    %
    %  Eye ball class for human, including human optics, sensor (cone
    %  mosaic, eye movement, macular) etc
    %
    %  Example:
    %    scene = sceneCreate;
    %    scene = sceneSet(scene, 'fov', 1);
    %    eb = eyeball(scene);
    %    visualize(eb, 'oi')
    %
    %    set(eb, 'sensor positions', randn(240, 2) * 5);
    %    compute(eb); % recompute cone absorptions
    %    visualize(eb, 'volts movie');
    %
    %
    %  HJ/BW, ISETBIO TEAM, 2015
    
    % private properties
    properties (GetAccess = public, SetAccess = private)
        oi
        sensor
    end
    
    % public methods
    methods (Access = public)
        % construction function
        % The inputs can be just parameters as name-value pairs or a scene
        % structure followed by name-value pairs
        %
        % Example:
        %   eb = eyeball();
        %   eb = eyeball('oi name', 'eyeball optical image');
        %   eb = eyeball(sceneCreate);
        %   eb = eyeball(sceneCreate, 'sensor name', 'eyeball sensor');
        function obj = eyeball(varargin)
            % check if scene is passed in
            if ~isempty(varargin) && isfield(varargin{1}, 'type') && ...
                    strcmp(varargin{1}.type, 'scene')
                % create default eyeball
                obj = eyeball(varargin{2:end});
                
                % adjust oi parameters and compute irradiance
                obj = compute(obj, varargin{1});
                
                return;
            end
            
            if mod(length(varargin), 2) ~= 0
                error('parameters should be name-value pairs');
            end
            
            obj.oi = oiCreate('human');
            obj.sensor = sensorCreate('human');
            
            for ii = 1 : 2 : length(varargin)
                set(obj, varargin{ii}, varargin{ii+1});
            end
        end
        
        % set function
        % In set funciton, we try to set to both oi and sensor
        % If just want to set to only oi or only sensor, format param as
        % 'oi param' or 'sensor param'
        %
        % Example:
        %   eb = eyeball();
        %   set(eb, 'oi name', 'human optical image');
        function obj = set(obj, param, val, varargin)
            [oType, p] = ieParameterOtype(param);
            if isequal(oType, 'sensor') % set to sensor
                if isempty(p)
                    if ~sensorCheckHuman(val), error('Invalid val'); end
                    obj.sensor = val;
                else
                    obj.sensor = sensorSet(obj.sensor, p, val,varargin{:});
                end
                
                
            elseif isequal(oType, 'oi') % set to oi
                if isempty(p)
                    if ~isequal(val.type, 'oi'), error('Invalid val'); end
                    obj.oi = val;
                else
                    obj.oi = oiSet(obj.oi, p, val, varargin{:});
                end
                
            else
                try
                    obj.oi = oiSet(obj.oi, p, val, varargin{:});
                    obj.sensor=sensorSet(obj.sensor,p,val,varargin{:});
                catch
                end
            end
        end
        
        % In get funciton, we try to get parameter values from eyeball
        % structure
        % If want to specify whether to get from oi or sensor,
        % format param as 'oi param' or 'sensor param'. Otherwise, we will
        % first try to get the parameter from OI and if failed, we try
        % again with sensor
        %
        % Example:
        %   eb = eyeball();
        %   set(eb, 'oi name', 'human optical image');
        function val = get(obj, param, varargin)
            [oType, p] = ieParameterOtype(param);
            if isequal(oType, 'sensor')
                if isempty(p)
                    val = obj.sensor;
                else
                    val = sensorGet(obj.sensor, p, varargin{:});
                end
            elseif isequal(oType, 'oi')
                if isempty(p)
                    val = obj.oi;
                else
                    val = oiGet(obj.oi, p, varargin{:});
                end
            else
                try
                    val = oiGet(obj.oi, p, varargin{:});
                catch
                    val = sensorGet(obj.sensor, p, varargin{:});
                end
            end
        end
        
        % Compute function
        % In this function, irradiance images and cone absorptions are
        % computed
        %
        % type can be strings choosen from 'oi', 'sensor' or 'both'. By
        % default both irradiance and absorption will be computed if scene
        % is supplied. If scene is not given, cone absorption get
        % recomputed
        function obj = compute(obj, scene, type)
            % check inputs
            if isempty(get(obj, 'oi photons'))
                if notDefined('scene'), error('scene required'); end
            elseif notDefined('scene')
                scene = [];
                if notDefined('type'), type = 'sensor'; end
            end
            if notDefined('type'), type = 'both'; end
            type = ieParamFormat(type);
            
            % compute irradiance
            if isequal(type, 'both') || isequal(type, 'oi')
                obj.oi = oiCompute(scene, obj.oi);
            end
            
            % adjust sensor parameters and compute absorptions
            if isequal(type, 'both')
                fov = sceneGet(scene, 'h fov');
                set(obj, 'sensor fov', fov, scene, obj.oi);
            end
            
            if isequal(type, 'both') || isequal(type, 'sensor')
                obj.sensor = coneAbsorptions(obj.sensor, obj.oi);
            end
        end
        
        % Visualize the computed data in eyeball object
        % We now can accept type as
        %   'oi'     - show OI window (default)
        %   'sensor' - show sensor window
        %   'volts movie'   - movie of volts sequences
        %   'adapted movie' - movie of adapted current
        %
        % For showing movie, additional parameter filename and framerate
        % can be specified
        function visualize(obj, type, varargin)
            % check inputs
            if notDefined('type'), type = 'oi'; end
            
            % visualize
            switch ieParamFormat(type)
                case 'oi'
                    vcAddObject(obj.oi); oiWindow(varargin{:});
                case 'sensor'
                    vcAddObject(obj.sensor); sensorWindow(varargin{:});
                case 'voltsmovie'
                    obj.showMovie(get(obj, 'sensor volts'), varargin{:});
                case 'adaptedmovie'
                    [~, adaptedCur] = coneAdapt(obj.sensor, 'rieke');
                    obj.showMovie(adaptedCur, varargin{:});
                otherwise
                    error('unknown type: %s', type);
            end
        end
        
    end
    
    methods (Access = private)
        % generate a movie sequence of the data
        % if fName is not empty, movie is written to file
        % if fName is empty, play the video directly
        function F = showMovie(~, data, fName, frameRate)
            % check inputs
            assert(ndims(data) == 3, 'data should be 3D');
            if notDefined('frameRate'), frameRate = 24; end
            nFrames = size(data, 3);
            
            % scale data
            data = ieScale(data, 0, 1);
            
            % play / write it to video file
            if notDefined('fName') % play the video
                implay(data, frameRate);
            else % write to file
                % generate movie
                vObj = VideoWriter(fName);
                vObj.FrameRate = frameRate;
                open(vObj);
                F = struct('cdata', [], 'colormap', gray(256));
                for ii = 1 : nFrames
                    F.cdata = uint8(255 * data(:,:,ii));
                    writeVideo(vObj, F);
                end
                close(vObj);
            end
        end
    end
end