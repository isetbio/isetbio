classdef cMosaic < handle
    
    properties (Constant)
        % Cone types
        LCONE_ID = 1;
        MCONE_ID = 2;
        SCONE_ID = 3;
        KCONE_ID = 4;
    end
    
    % Public properties
    properties
        eccentricityDegs;
        sizeDegs;
        whichEye;
        coneDensities;
        tritanopicRadiusDegs;
        noiseFlag;
        
        % Color for cone rendering in mosaic visualization
        lConeColor = [1 0.5 0.6];
        mConeColor = [0.2 1 0.5];
        sConeColor = [0.5 0.1 1];
        kConeColor = [0.9 0.9 0.2];
    end
    
    % Read-only properties
    properties (GetAccess=public, SetAccess=private)
        % [n x 2] matrix of RGC positions, in microns
        coneRFpositionsMicrons;
        
        % [n x 2] matrix of RGC positions, in degrees
        coneRFpositionsDegs;
        
        % [n x 1] vector of cone spacings
        coneRFspacingsDegs;
        coneRFspacingsMicrons;
        
        % [n x 1] vector of cone types
        coneTypes;
        
        % indices of cones
        lConeIndices;
        mConeIndices;
        sConeIndices;
        kConeIndices;
    end
    
    % Private properties
    properties (GetAccess=private, SetAccess=private)
        % Size (in degs) of source lattice from which to crop positions for
        % the desired eccentricity
        sourceLatticeSizeDegs = 45;
    end
    
    methods
        function obj = cMosaic(varargin)
            
            % Parse input
            p = inputParser;
            p.addParameter('eccentricityDegs', [0 0], @(x)(isnumeric(x) && (numel(x) == 2)));
            p.addParameter('sizeDegs', [0.2 0.2], @(x)(isnumeric(x) && (numel(x) == 2)));
            p.addParameter('whichEye', 'right', @(x)(ischar(x) && (ismember(x, {'left', 'right'}))));
            p.addParameter('coneDensities', [0.6 0.3 0.1 0.0], @isnumeric);
            p.addParameter('tritanopicRadiusDegs', 0.15, @isscalar);
            p.addParameter('noiseFlag', 'random', @(x)(ischar(x) && (ismember(x, {'random', 'frozen', 'none'}))));
            p.parse(varargin{:});
            
            obj.eccentricityDegs = p.Results.eccentricityDegs;
            obj.sizeDegs = p.Results.sizeDegs;
            obj.whichEye = p.Results.whichEye;
            obj.coneDensities = p.Results.coneDensities;
            obj.tritanopicRadiusDegs = p.Results.tritanopicRadiusDegs;
            obj.noiseFlag = p.Results;
            
            % Initialize positions
            obj.initializeConePositions();
            
            % Assign types
            obj.assignConeTypes();
        end
        
        % Method to visualize the cone mosaic
        visualize(obj, varargin);
    end
    
    methods (Access=private)
        initializeConePositions(obj);
        assignConeTypes(obj);
    end
    
    methods (Static)
        % Method to render an array of patches
        renderPatchArray(axesHandle, pixelOutline, xCoords, yCoords, faceColors, edgeColor, lineWidth);
    end
    
end

