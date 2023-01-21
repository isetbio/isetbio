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

        rgcsNum;
        inputConesNum;

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

            % Parse input
            p.parse(sourceMidgetRGCMosaic, varargin{:});
            

            % The sourceMidgetRGCMosaic must be frozen
            assert(sourceMidgetRGCMosaic.isFrozen, 'SourceMidgetRGCMosaic must be frozen');

            obj.name = p.Results.name;
            obj.sizeDegs = p.Results.sizeDegs;
            
            if (~isempty(p.Results.eccentricityDegs))
                % eccentricityDegs is always used if present
                obj.eccentricityDegs = p.Results.eccentricityDegs;
            elseif (~isempty(p.Results.positionDegs))
                % Use positionDegs if available but there is not
                % eccentricityDegs
                obj.eccentricityDegs = p.Results.positionDegs;
            end

            % Generate the mRGCMosaic by cropping the sourceMidgetRGCMosaic
            obj.generateByCroppingTheSourceMosaic(p.Results.sourceMidgetRGCMosaic);

        end % Constructor

        % Compute method
        [response, responseTemporalSupport] = compute(obj, ...
            theConeMosaicResponse, theConeMosaicResponseTemporalSupport, varargin);

    end % Public methods

    % Private methods
    methods (Access=private)
        % Method to generate the mRGCMosaic by cropping the sourceMidgetRGCMosaic
        generateByCroppingTheSourceMosaic(obj, sourceMidgetRGCMosaic);
    end % Private methods

    
end
