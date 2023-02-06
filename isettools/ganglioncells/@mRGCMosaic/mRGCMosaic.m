classdef mRGCMosaic < handle
    % Create an mRGCMosaic for computation by cropping a subregion from a
    % previously frozen midgetRGCMosaic (which takes hours/days to generate)
    %
    %
    % Syntax:
    %   theMRGCMosaic = mRGCMosaic(sourceMidgetRGCMosaic);
    %   theConeMosaic = cMosaic(sourceMidgetRGCMosaic, 'sizeDegs',[5,5],'eccentricityDegs',[0,0]);
    %
    % Description:
    %    Creates 
    %
    % Inputs:
    %    None required.
    %
    % Outputs:
    %    mRGCMosaic                 The created mRGCMosaic object.
    %
    % Optional key/value pairs:
    %    'name'                     - String. Mosaic name. 
    %    'eccentricityDegs'         - [x,y]. Center of the mosaic, in degrees. Default: [0 0]
    %    'sizeDegs'                 - [sx, sy]. Width and height of the mosaic, in degrees. Default: [0.5 0.5]
    %
    %
    %
    % History:
    %    January 2023  NPC  Wrote it

    % Public properties
    properties  (GetAccess=public, SetAccess=public)
        % Name of the mosaic
        name;

        % Boolean. Set to true to allow the object to print various info
        beVerbose;
    end

    % Read-only properties
    properties (GetAccess=public, SetAccess=private)
        
        % Eccentricity [x,y] of the center of the mRGC mosaic (degs)
        eccentricityDegs;

        % Size of the mRGC mosaic (degs)
        sizeDegs;

        % The input cone mosaic (@cMosaic)
        inputConeMosaic;

        % [n x 2] matrix of RF positions, in degrees
        rgcRFpositionsDegs;

        % The (sparse) [conesNum x rgcRFsNum] conePoolingMatrices for the
        % center and the surround subregions
        centerConePoolingMatrix;
        surroundConePoolingMatrix;

        % The type of the majority of cone inputs to the RF center subregions
        centerSubregionMajorityConeTypes;

        % How many center cones are of type cMosaic.LCONE_ID, cMosaic.MCONE_ID, cMosaic.SCONE_ID
        centerSubregionConeTypesNums;

        % How many cones in each RF center
        centerSubregionConesNums;

        % The number of RGCs
        rgcsNum;

        % The number of input cones
        inputConesNum;

        % The spatial position grid at which the multifocal RTVF model was fitted
        multifocalRTVFgrids;

        % The params for generating optics at the positions where the
        % multifocal RTVF model was fitted
        multifocalRTVFopticsParams;
    end

    % Private properties
    properties (GetAccess=private, SetAccess=private)
        % Cache for various visualizations
        visualizationCache = [];
    end

    % Public methods
    methods
        
        % Constructor
        function obj = mRGCMosaic(sourceMidgetRGCMosaic, varargin)
            % Parse input
            p = inputParser;
            p.addRequired('sourceMidgetRGCMosaic', @(x)(isa(x, 'midgetRGCMosaic')));
            p.addParameter('name', 'my midget RGC mosaic', @ischar);
            p.addParameter('eccentricityDegs', [], @(x)(isempty(x))||(isnumeric(x) && (numel(x) == 2)));
            p.addParameter('positionDegs', [], @(x)(isempty(x))||(isnumeric(x) && (numel(x) == 2)));
            p.addParameter('sizeDegs', [], @(x)(isempty(x))||(isnumeric(x) && (numel(x) == 2)));
            p.addParameter('visualizeSpatialRelationshipToSourceMosaic', false, @islogical);
            p.addParameter('beVerbose', false, @islogical);
            % Parse input
            p.parse(sourceMidgetRGCMosaic, varargin{:});
            

            % The sourceMidgetRGCMosaic must be frozen
            assert(sourceMidgetRGCMosaic.isFrozen, 'sourceMidgetRGCMosaic must be frozen.');

            obj.name = p.Results.name;
            obj.sizeDegs = p.Results.sizeDegs;
            obj.beVerbose = p.Results.beVerbose;
            
            if (~isempty(p.Results.eccentricityDegs))
                % eccentricityDegs is always used if present
                obj.eccentricityDegs = p.Results.eccentricityDegs;
            elseif (~isempty(p.Results.positionDegs))
                % Use positionDegs if available but there is not
                % eccentricityDegs
                obj.eccentricityDegs = p.Results.positionDegs;
            end

            % Get the input cone mosaic
            obj.inputConeMosaic = sourceMidgetRGCMosaic.inputConeMosaic;

            % Get the optics params used to train the multi-focal RTVF
            % so we can pass them to the user when requested
            obj.retrieveMultifocalRTVFOpticsParams(...
                p.Results.sourceMidgetRGCMosaic);

            % Generate the mRGCMosaic by cropping the sourceMidgetRGCMosaic
            obj.generateByCroppingTheSourceMosaic(...
                p.Results.sourceMidgetRGCMosaic, ...
                p.Results.visualizeSpatialRelationshipToSourceMosaic);

            % The type of the majority of cone types in the RF center
            [~, ...
                obj.centerSubregionConeTypesNums, ...
                obj.centerSubregionMajorityConeTypes, ...
                obj.centerSubregionConesNums] = obj.centerConeTypeWeights(1:obj.rgcsNum);
            
        end % Constructor

        % Method to compute optics at a given position
        theOI = multiFocalRTVFopticsAtPosition(obj, eccDegs);

        % Compute method
        [response, responseTemporalSupport] = compute(obj, ...
            theConeMosaicResponse, theConeMosaicResponseTemporalSupport, varargin);

        % Method to visualize the mRGC mosaic and its activation (RF centers only)
        visualize(obj, varargin);

        % Method to visual the RF structure of individual RGCs
        visualizeRFs(obj, rgcIndices, varargin);

        % Method to return the net weights and numbers of types of cones in the RF center
        % as well as the type with the strongest net weight. If 2 types
        % of cones have net weights that are close to each other (90%)
        % then theMajorityConeType is set to nan
        [theCenterConeTypesWeights, theCenterConeTypesNum, theMajorityConeTypes, theCenterConesNum] = ...
            centerConeTypeWeights(obj, theRGCindices);

    end % Public methods

    % Private methods
    methods (Access=private)
        % Method to generate the mRGCMosaic by cropping the sourceMidgetRGCMosaic
        generateByCroppingTheSourceMosaic(obj, ...
            sourceMidgetRGCMosaic, ...
            visualizeSpatialRelationshipToSourceMosaic);

        % Method to retrieve the optics params of the multifocal RTVFobject
        % used to derive cone pooling weights
        retrieveMultifocalRTVFOpticsParams(obj, sourceMidgetRGCMosaic);
    end % Private methods

    
end
