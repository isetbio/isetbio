classdef rgcParameters < handle
%% Creates a class definition (structure) of the RGC parameters
%
%   rgcParameters
%
% The RGC parameters creation makes most of what you need for an RGC
% simulation, but it does not create a layer.  Thus, you must follow the
% parameter creation call with a call to add layers. See the example below.
%
% To see the list of rgcParameters, use:
%
%   properties('rgcParameters')   or methodsview('rgcParameters')
%
% When you edit this file, you must clear the class definitions before
% re-running functions or scripts that use this class.  To reset, use
%
%   clear rgcParameters
%
% More explanations of the  rgcParameters are below. The functions that
% implement getting and setting most of the parameters are rgcGet and
% rgcSet.m. Some of the functions are also implemented here.
%
% Discussion
%
%  Talk about how rgcP contains scene, oi, and sensor data Talk about how
%  to set and get parameters Talk about the compute routines
%
% Examples:
%
%   rgcP = rgcParameters;
%
% You can create several type of layers, using rgcP.addLayer(cellType)
% where the parameters for cellType are established in layerSetCellType.
% Current options include 'on midget','off midget', and so forth.
%
%  rgcP.addLayer('on midget');
%
% See the examples in s_rgcScene2Cones and s_rgcCones2RGC for examples of
% how to set the properties.
%
%
% (c) 2010 Stanford Synapse Team
    
%% Initialize the rgcParameters object
%
% We are organized in terms of public, private, methods
% I am not sure why the methods are broken up instead of grouped into one
% big methods block.
% Also, note that there is a function at the bottom.
%
    
    % These parameters can be set or get by any function
    properties(GetAccess = 'public', SetAccess = 'public')
        noise    = defaultNoise();
        
        name = sprintf('rgcSimulation-%s',date);  % Simulation name
        
        absorptions = [];  % 3D array of cone absorptions (x,y,t)
        cVolts = [];       % Cone voltages (x,y,t). Should replace absorptions.
        scene  = [];       % ISET scene
        oi     = [];       % Optical image
        
        % We always need a sensor to determine some parameters. Typically,
        % this field is over-written by the calling script.
        sensor = sensorCreate('human');
        
        % The temporal response is slow, so we always run to about 150 ms
        trDur = 150; % Input response function duration (ms)
        
        % Change these parameters to be named noiseMeanV and noiseSDV Or, just
        % get rid of them.  We should add rgc noise using a separate noise
        % process
        %
        % meanV = 0.002; % The mean of the RGC noise (volts)
        % stdV  = 0.001; % The standard-deviation of the RGC noise (volts)
        layers      =  {};   % Empty cell array of layers
    end
    
    % Private get and set properties
    % These can only be set and get by the rgcP object itself.
    properties(GetAccess = 'private', SetAccess = 'private')
        % layer properties, default means no layer
        
    end
    
    % These can be read by any routine, but only set by this object
    properties(GetAccess = 'public', SetAccess = 'private')
        
    end
    
    % These parameters can be derived from the parameters stored in the
    % rgcP object.  These parameters are not stored directly.
    properties(Dependent = true, SetAccess = 'private')
        %     dT;                 % sampling time step between frames (ms)
        %     coneSpacing;        % dS between cones in um.
        %    feedbackTimes;      % sampled times to apply to feedback response
        %    coneGrid;           % Locations of the cones
        %     nLayers;            % number of layers
        distanceFunction;   % distance function for the coupling
    end
    
    
    
    %% Set methods
    % We work through rgc"Get/Set"Parameters calls that include aliases for parameters
    % and management of the string format.
    methods
        % General set parameters function
        function res = set(obj,paramName,value,varargin)
            res = rgcSet(obj,paramName,value,varargin);
        end
        
        % Set absorptions - either from a structure or as a matrix
        function set.absorptions(obj,value)
            % The absorptions were initially used to store cone voltages.  We
            % are trying to separate these out now.  In the future, the
            % absorptions may not be stored here, only the voltages
            % Other times, we just read a 3D array of absorptions
            obj.absorptions = value;
            
        end
        
        function set.cVolts(obj,value)
            % The absorptions were initially used to store cone voltages.  We
            % are trying to separate these out now.  In the future,
            % we may only include the voltages
            % Other times, we just read a 3D array of absorptions
            obj.cVolts = value;
        end
        
        % Set ISET sensor
        function set.sensor(obj,value)
            obj.sensor   = value;
        end
        
        % Set ISET optical image structure
        function set.oi(obj,value)
            obj.oi   = value;
        end
        
    end
    
    %% Get methods
    % Mainly we work through rgc"Get/Set"Parameters calls that include aliases
    % for parameters and management of the string format.
    methods
        % This is the class get function.  We have a couple of special cases
        % (below), but in the future there should probably be only the one get
        % function, managed here through rgcGet().
        function res = get(obj,paramName,optParam)
            if notDefined('optParam'), res = rgcGet(obj,paramName);
            else res = rgcGet(obj,paramName,optParam);
            end
        end
        
    end
    %% Additional get methods for dependent parameters
    %
    methods
        function distfunc = get.distanceFunction(obj)
            distfunc = @(D, wRange) wRange * exp(-D );
        end
        
    end
    
    %% Methods for rgcP displaying, copying, resetting, cleaning...
    methods
        
        % Constructor for the rgcParameters class
        function obj = rgcParameters
            obj.sensor = sensorSet(obj.sensor,'exp time',0.001);
        end
        
        % Display parameters
        function Display(obj,paramName,varargin)
            % rgcP.Display;
            if notDefined('paramName'), paramName = 'name'; end
            if notDefined('varargin'), rgcDisplayParameter(obj,paramName);
            else                       rgcDisplayParameter(obj,paramName,varargin);
            end
        end
        
        %
        function newObj = Copy(obj)
            % r = rgcP.Copy;
            newObj = rgcParameters;
            %public values
            rgcCopy(obj,newObj); % copies obj public values in newObj values
            %layers
            nL = obj.get('nLayers');
            layersC = cell(1,nL);
            for ii = 1:nL
                layersC{ii} = obj.getLayer(ii).Copy();
                layersC{ii}.parent = newObj;
            end
            newObj.layers = layersC;
            % private values
            newObj.sensor = obj.sensor;
            newObj.oi = obj.oi;
        end
        
        % Reset object defaults
        function newObj = resetDefaultsParameters(obj)
            % r = rgcP.resetDefaultsParameters;
            newObj = rgcParameters();
            %resetting the layers
            nL = obj.get('nLayers');
            for i = 1:nL
                newObj = newObj.addLayer();
            end
            % Preserve the absorptions
            newObj.cVolts = obj.cVolts;
            rgcCopy(newObj,obj); % copies newObj values in obj values
        end
        
    end
    
    %% Methods for managing layers
    methods
        % Retrieve a layer
        function layer = getLayer(obj,n)
            nL = obj.get('nLayers');
            if notDefined('n')
                if ( nL == 0), error('One layer required: Use rgcP.addLayer();');
                elseif ( nL == 1), n = 1;
                else    error('Specify a layer, e.g., param.getLayer(2);');
                end
            end
            if (n<1 || n>nL), error('getLayer(obj,n): n has to be in [|0 %d|]',nL); end
            layer = obj.layers{n};
        end
        
        % Adding one layer
        function obj = addLayer(obj,cellType,ang)
            % Add a layer of ganglion cells to the rgcP structure
            % The properties of the different types are stored in the function
            % layerSetCellType.
            if notDefined('cellType')
                cellType = 'default';
                disp('Using "default" layer type.');
            end
            if notDefined('ang'), ang = []; end
            
            % Creating the layer
            layer = rgcLayer(obj);
            
            % Setting it to the wanted type
            layerSetCellType(layer,cellType,ang);
            
            % putting it into the layers field
            obj.layers{obj.get('nLayers') + 1} = layer;
        end
        
        %Removing a layer
        function obj = removeLayers(obj,nList)
            nL = obj.get('nLayers');
            if (min(nList)<1 || max(nList)>nL)
                error('The layer numbers (%.0f) are incorrect',nList);
            end
            obj.layers(nList)=[];
        end
        
        %Switching two layers
        function obj = switchLayers(obj,i,j)
            nL = obj.get('nLayers');
            if (i<1 || i>nL || j<1 || j>nL)
                error('The layer numbers (%.0f, %.0f) are incorrect',i,j);
            end
            tmpLayer = obj.layers{i};
            obj.layers{i} = obj.layers{j};
            obj.layers{j} = tmpLayer;
        end
    end
    
end  % rgcP methods.

%% Functions for creating default structures
function noise = defaultNoise()
% Noise function and function parameters
% Gaussian noise with size S, mean m and standard deviation s
noise.nsFunc = @(S, s, m)randn(S) * s + m;
noise.nFrame = 0;        % Number of different noise frames
end