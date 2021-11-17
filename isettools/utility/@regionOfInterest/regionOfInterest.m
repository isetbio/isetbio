classdef regionOfInterest < handle
    % Create an regionOfInterest object
    %
    % Syntax:
    %   rectangularROI = regionOfInterest(...
    %       struct(...
    %           'units', 'degs', ...
    %           'shape', 'rect', ...
    %           'center', [14, 2], ...
    %           'width', 10, ...
    %           'height', 5, ...
    %           'rotation', 30.0...
    %       ));
    %
    %   ellipticalROI = regionOfInterest(...
    %       struct(...
    %           'units', 'degs', ...
    %           'shape', 'ellipse', ...
    %           'center', [14, 2], ...
    %           'majorAxisDiameter', 10, ...
    %           'minorAxisDiameter', 5, ...
    %           'rotation', 30.0...
    %       ));
    %
    % Description:
    %    A @regionOfInterest object handles the geometry of an ROI and
    %    provides various convenience functions for visualizing the ROI, 
    %    for detecting points inside/outside the ROI, etc.
    %
    %
    % Inputs:
    %    None required.
    %
    % Outputs:
    %    regionOfInterest           The created regionOfInterest object.
    %
    % See:
    %   t_roiBasic
    %

    % History:
    %    November 2021  NPC  Wrote it
    
    % Read-only properties
    properties (GetAccess=public, SetAccess=private)
        
        % Struct defining the geometry of the regionOfInterest
        geometryStruct
    end
    
    properties (Dependent)
        outline
    end
    
    % Public methods
    methods
        % Constructor
        function obj = regionOfInterest(geometryStruct, varargin)
             % parse input
            p = inputParser;
            p.addRequired('geometryStruct', @isstruct);
            p.parse(geometryStruct, varargin{:});
            
            obj.validateGeometry(p.Results.geometryStruct);
        end
        
        % Method to visualize the ROI
        visualize(obj, varargin);
        
        % Method to return the indices of the points that lie within the ROI
        idx = indicesOfPointsInside(obj, points)
        
        % Method to return the indices of the points that lie outside of the ROI
        idx = indicesOfPointsOutside(obj, points)
        
        % Getter for outline
        function val = get.outline(obj)
            val = generateOutline(obj);
        end
    end
    
    
    methods (Static = true)
    end
    
    methods (Access=private)
        validateGeometry(obj, geoStruct);
    end
    
end
