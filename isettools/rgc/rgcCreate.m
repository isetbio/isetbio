classdef rgcCreate < handle
%% Creates a class definition of retina ganglion cell parameters
%
%   rgc = rgcCreate;
%
% The rgcCreate makes most of what you need for an RGC simulation, but it
% does not create a layer.  Thus, you must follow the parameter creation
% call with a call to add layers.
%
% To see the list of rgc properties and methods, use:
%
%   properties('rgcCreate')   or methodsview('rgcCreate')
%
% Notes:
%   This file is adopted from rgcParameters.m, written by BW / JW
%
% Examples:
%
%   rgcP = rgcParameters;
%
% You can create several type of layers, using rgcP.addLayer(cellType)
% where the parameters for cellType are established in layerSetCellType.
% Current options include 'on midget', 'off midget', and so forth.
%
%  rgcP.addLayer('on midget');
%
% See also:
%   rgcSet, rgcGet
%
% (HJ) ISETBIO TEAM, 2015
    
%% Initialize the rgcParameters object

    %  public properties
    properties(GetAccess = 'public', SetAccess = 'public')
        name = sprintf('rgc-%s',date);  % Simulation name
        
        % Charaterize noise
        % Gaussian noise with size S, mean m and standard deviation s
        noise = @(S, s, m)randn(S) * s + m;
        
        % The temporal response is slow, so we always run to about 150 ms
        trDur  = 150;   % Input response function duration (ms)
        layers =  {};   % Empty cell array of layers
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
        % dT;               % sampling time step between frames (ms)
        % coneSpacing;      % dS between cones in um.
        % feedbackTimes;    % sampled times to apply to feedback response
        % coneGrid;         % Locations of the cones
        % nLayers;          % number of layers
        distanceFunction;   % distance function for the coupling
    end
    
    
    
%% Get / Set methods
% We work through rgcGet/ rgcSet calls that include aliases for parameters
% and management of the string format
    methods
        % This is the class set function, see rgcSet for details
        function res = set(obj,paramName,value,varargin)
            res = rgcSet(obj,paramName,value,varargin);
        end
        
        % This is the class get function, see rgcGet for details
        function res = get(obj,paramName,optParam) 
            if notDefined('optParam'), res = rgcGet(obj,paramName);
            else res = rgcGet(obj,paramName,optParam);
            end
        end
    end
    
%% Methods for rgc displaying, copying, resetting and cleaning
    methods
        % Display parameters
        function Display(obj,paramName,varargin)
            % rgc.Display;
            if notDefined('paramName'), paramName = 'name'; end
            if notDefined('varargin'), rgcDisplayParameter(obj,paramName);
            else rgcDisplayParameter(obj,paramName,varargin);
            end
        end
        
        %
        function newObj = Copy(obj)
            % r = rgcP.Copy;
            newObj = rgcCreate;
            %public values
            rgcCopy(obj,newObj); % copies public values in newObj values
            %layers
            nL = obj.get('nLayers');
            layersC = cell(1,nL);
            for ii = 1:nL
                layersC{ii} = obj.getLayer(ii).Copy();
                layersC{ii}.parent = newObj;
            end
            newObj.layers = layersC;
        end
    end
    
%% Methods for managing layers
    methods
        % Retrieve a layer
        function layer = getLayer(obj,n)
            nL = obj.get('nLayers');
            if notDefined('n')
                if nL == 0
                    error('One layer required: Use rgcP.addLayer();');
                elseif nL == 1, n = 1;
                else   error('Specify a layer, e.g., param.getLayer(2);');
                end
            end
            if n<1 || n>nL
                error('getLayer(obj,n): n has to be in [|0 %d|]',nL);
            end
            layer = obj.layers{n};
        end
        
        % Adding one layer
        function obj = addLayer(obj, cellType, ang)
            % Add a layer of ganglion cells to the rgcP structure
            % The properties of the different types are stored in the 
            % function layerSetCellType
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
        function obj = removeLayers(obj, nList)
            nL = obj.get('nLayers');
            if min(nList)<1 || max(nList)>nL
                error('The layer numbers (%.0f) are incorrect', nList);
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
end