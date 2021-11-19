classdef regionOfInterest < handle
    % Create an regionOfInterest object
    %
    % Syntax:
    %   rectangularROI = regionOfInterest('geometryStruct', ...
    %       struct(...
    %           'units', 'degs', ...
    %           'shape', 'rect', ...
    %           'center', [14, 2], ...
    %           'width', 10, ...
    %           'height', 5, ...
    %           'rotation', 30.0...
    %       ));
    %
    %    regionOfInterest('help');
    %    regionOfInterest('params');
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
    
    % Constant properties
    properties (Constant)
        validShapes = {'rect', 'ellipse', 'line'};
        validUnits = {'degs', 'microns'};
    end
    
    % Private properties
    properties (GetAccess=private, SetAccess=private)   
        defaultGeometryStruct = struct(...
                    'units', 'degs', ...
                    'shape', 'rect', ...
                    'center', [0, 0], ...
                    'width', 1, ...
                    'height', 1, ...
                    'minorAxisDiameter', 1, ...
                    'majorAxisDiameter', 1, ...
                    'from', [-1 0], ...
                    'to', [1 0], ...
                    'rotation', 0.0...
                    );     
    end
    
    % Public properties
    properties (Access=public)
        units
        shape
        center
        width
        height
        minorAxisDiameter
        majorAxisDiameter
        from
        to
        rotation
    end
    
    % Dependent 
    properties (Dependent)
        outline
        geometryStruct
    end
    
    % Public methods
    methods
        % Constructor
        function obj = regionOfInterest(varargin)
            % parse input
         
            if (nargin == 1)
                mode = varargin{1};
                switch (mode)
                    case 'help'
                        regionOfInterest.displayHelp();
                    case 'params'
                        visualizeParamsStructTree(obj.defaultGeometryStruct, 'default params');
                    otherwise
                        fprintf(2,'Incorrect instantiation. See below for ways to instantiate a @regionOfInterest.\n');
                        regionOfInterest.displayHelp();
                        error('Incorrect instantiation');
                end
                return;
            end
           
            % Parse input
            p = inputParser;
            p.addParameter('units',    obj.defaultGeometryStruct.units,  @(x)(ischar(x) && (ismember(x, regionOfInterest.validUnits))));
            p.addParameter('rotation', obj.defaultGeometryStruct.rotation, @isscalar);
            p.addParameter('shape',    obj.defaultGeometryStruct.shape, @(x)(ischar(x) && (ismember(x, regionOfInterest.validShapes))));
            p.addParameter('center',   obj.defaultGeometryStruct.center, @(x)(isnumeric(x) && ((isempty(x))||(numel(x) == 2))));
            p.addParameter('width',    obj.defaultGeometryStruct.width, @isscalar);
            p.addParameter('height',   obj.defaultGeometryStruct.height, @isscalar);
            p.addParameter('minorAxisDiameter', obj.defaultGeometryStruct.minorAxisDiameter, @isscalar);
            p.addParameter('majorAxisDiameter', obj.defaultGeometryStruct.majorAxisDiameter, @isscalar);
            p.addParameter('from',    obj.defaultGeometryStruct.from, @isscalar);
            p.addParameter('to',      obj.defaultGeometryStruct.to, @isscalar);
            p.addParameter('geometryStruct', []);
            p.parse(varargin{:});
            
            if (~isempty(p.Results.geometryStruct))
                % Validate the input geometry struct
                obj.validateGeometry(p.Results.geometryStruct);
            else
            
                % Update default values based on key/value pairs
                fNames = setdiff(fieldnames(p.Results), 'geometryStruct');
                for k = 1:numel(fNames)
                    obj.(fNames{k}) = p.Results.(fNames{k});
                end
            end
        end
        
        % Global setter
        function set(obj, varargin)
            settableParams = fieldnames(obj.defaultGeometryStruct);
            settableParams{numel(settableParams)+1} = 'geometryStruct';
            for k = 1:2:numel(varargin)-1
                if (~ismember(varargin{k}, settableParams))
                    fprintf(2,'Field ''%s'' does not exis and is therefore not settable. Ignoring it.\n', varargin{k});
                else
                    obj.(varargin{k}) = varargin{k+1};
                end
            end
        end
        
        % Setters
        function set.units(obj,val)
            assert(ismember(val, regionOfInterest.validUnits));
            obj.units = val;
        end
        
        function set.shape(obj,val)
            assert(ismember(val, regionOfInterest.validShapes));
            switch (val)
                case 'line'
                    % Ensure that 'from' and 'to' are not-empty, and if
                    % they are set it to their default values
                    if (isempty(obj.from))
                        obj.from = obj.defaultGeometryStruct.from;
                    end
                    
                    if (isempty(obj.to))
                        obj.to = obj.defaultGeometryStruct.to;
                    end

                case 'ellipse'
                    % Ansure that 'rotation', 'minorAxisDiameter' and 'majorAxisDiameter' are not-empty, and if
                    % they are set it to their default values
                    if (isempty(obj.minorAxisDiameter))
                        obj.minorAxisDiameter = obj.defaultGeometryStruct.minorAxisDiameter;
                    end
                    if (isempty(obj.majorAxisDiameter))
                        obj.majorAxisDiameter = obj.defaultGeometryStruct.majorAxisDiameter;
                    end
                    if (isempty(obj.rotation))
                        obj.rotation = obj.defaultGeometryStruct.rotation;
                    end
                case 'rect'
                    % Ansure that 'rotation', 'width' and 'height' are not-empty, and if
                    % they are set it to their default values
                    if (isempty(obj.width))
                        obj.width = obj.defaultGeometryStruct.width;
                    end
                    if (isempty(obj.height))
                        obj.height = obj.defaultGeometryStruct.height;
                    end
                    if (isempty(obj.rotation))
                        obj.rotation = obj.defaultGeometryStruct.rotation;
                    end
                    
            end
            obj.shape = val;
        end
        
        % Setter for dependent geometryStruct
        function set.geometryStruct(obj,val)
            fNames = fieldnames(val);
            for k = 1:numel(fNames)
                obj.(fNames{k}) = val.(fNames{k});
            end
        end
        
        
        % Getter for dependent geometryStruct
        function val = get.geometryStruct(obj)
            val = struct(...
                'units', obj.units, ...
                'shape', obj.shape, ...
                'center', obj.center, ...
                'width', obj.width, ...
                'height', obj.height, ...
                'minorAxisDiameter', obj.minorAxisDiameter, ...
                'majorAxisDiameter', obj.majorAxisDiameter, ...
                'from', obj.from, ...
                'to', obj.to, ...
                'rotation', obj.rotation...
                );
        end
        
        % Getter for outline
        function val = get.outline(obj)
            val = generateOutline(obj);
        end
        
        
        % Method to visualize the ROI
        visualize(obj, varargin);
        
        % Method to return the indices of the points that lie within the ROI
        idx = indicesOfPointsInside(obj, points)
        
        % Method to return the indices of the points that lie outside of the ROI
        idx = indicesOfPointsOutside(obj, points)
        
        % Method to return the indices of the points that lie up to maxDistance from
        % the perimeter of the ROI, with the perimeter outline sampled using nSamples
        idx = indicesOfPointsAround(obj, points, pointsPerSmaple, samplingPoints, maxDistance);
    end
    
    
    methods (Static = true)
        function displayHelp()
            fprintf('---------------------------');
            fprintf('\n@regionOfInterest: A class to manage a region of interest.');
            fprintf('\nUsage patterns:');
            
            fprintf(2,'\n\n0. Print some documentation:');
            fprintf('\n\t regionOfInterest(''help'')');
            
            fprintf(2,'\n\n1. Generate the default ROI and visualize it:');
            fprintf('\n\t d = regionOfInterest(); d.visualize();');
            
            fprintf(2,'\n\n2. Generate default ROI with some modifications and visualize it:');
            fprintf('\n\t d = regionOfInterest(''center'', [3 5], ''width'', 3, ''rotation'', 30); d.visualize();');
            
            fprintf(2,'\n\n3. Alter some params of an existing ROI and visualize it:');
            fprintf('\n\t d.set(''center'', [-3 -5], ''height'', 5); d.visualize();');
            
            fprintf(2,'\n\n4. Generate a completely custom ROI and visualize it:');
            fprintf('\n\t d.set(''geometryStruct'', struct(...');
            fprintf('\n\t\t ''units'', ''degs'', ...');
            fprintf('\n\t\t ''shape'', ''ellipse'', ...');
            fprintf('\n\t\t ''center'', [15.5, 1.5], ...');
            fprintf('\n\t\t ''minorAxisDiameter'', 10, ...');
            fprintf('\n\t\t ''majorAxisDiameter'', 5, ...');
            fprintf('\n\t\t ''rotation'', 13.0));');
            fprintf('\n\t d.visualize();');

            fprintf(2,'\n\n5. Compute the indices of points that lie inside an ROI:');
            fprintf('\n\t randomPoints = bsxfun(@plus, [14 2], randn(600,2)*3);');
            fprintf('\n\t idxPointsInside = d.indicesOfPointsInside(randomPoints);');
            
            fprintf(2,'\n\n6. Compute the indices of points that lie outside an ROI:');
            fprintf('\n\t idxPointsOutside = d.indicesOfPointsOutside(randomPoints);');
            
            fprintf(2,'\n\n7. Compute the indices of points that lie near the perimeter of an ROI:');
            fprintf('\n\t samplingPoints = 1000; %% sample the perimeter of the ROI along 1000 points');
            fprintf('\n\t pointsPerSample = 30;  %% return up to 30 points for each sample along the perimeter');
            fprintf('\n\t maxDistance = 0.5;     %% points must be no further than 0.5 units away from the closest perimeter sample');
            fprintf('\n\t idx = d.indicesOfPointsAround(randomPoints, pointsPerSample, samplingPoints, maxDistance);');
            fprintf('\n---------------------------\n');
        end
    end
    
    methods (Access=private)
        validateGeometry(obj, geoStruct);
    end
    
end
