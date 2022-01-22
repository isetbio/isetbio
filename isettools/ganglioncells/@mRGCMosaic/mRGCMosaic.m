classdef mRGCMosaic < handle
    % Create an midget RGC mosaic object
    %
    % Syntax:
    %   theMRGCmosaic = mRGCMosaic;
    %
    % Description
    %
    %
    % Inputs:
    %    None required.
    %
    % Outputs:
    %    mRGCMosaic           The created mRGCMosaic object.
    %
    % Optional key/value pairs:
    %    'name'                                 - String. Mosaic name. Default: 'mRGCmosaic'.
    %    'centerPolarity'                       - String. Center polarity.  Valid options are: {'ON', 'OFF'}. Default: 'ON'.
    %    'whichEye'                             - Which eye. Valid options are: {'left eye', 'right eye'}. Default: 'right eye'.
    %    'eccentricityDegs'                     - [x,y]. Center of the mosaic, in degrees. Default: [0 0]
    %    'sizeDegs'                             - [sx, sy]. Width and height of the mosaic, in degrees. Default: [0.4 0.4]
    %    'computeMeshFromScratch'               - Logical, indicating whether to recompute the mesh from scratch. Default: false
    %    'visualizeMeshConvergence'             - Logical, indicating whether to visualize the convergence of the mesh. Default: false 
    %    'exportMeshConvergenceHistoryToFile'   - Logical, indicating whether to save the convergence of the mesh. Default:false
    %    'maxMeshIterations'                    - Scalar. Max number of iterations for the mesh generation. Default: 100.
    
    % Public properties
    properties  (GetAccess=public, SetAccess=public)
        % Name of the mosaic
        name;

        % Noise flag for RGCs
        noiseFlag;
        
        % Random seed
        randomSeed;
    end

    % Read-only properties
    properties (GetAccess=public, SetAccess=private)

        % [n x 2] matrix of mRGC RF positions, in microns
        rgcRFpositionsMicrons;
        
        % [n x 2] matrix of mRGC RF positions, in degrees
        rgcRFpositionsDegs;
        
        % [n x 1] vector of mRGC RF spacings, in degrees
        rgcRFspacingsDegs;
        
        % [n x 1] vector of mRGC RF spacings, in microns
        rgcRFspacingsMicrons;

        % Eccentricity [x,y] of the center of the mRGC mosaic (degs)
        eccentricityDegs;

        % Eccentricity [x,y] of the center of the mRGC mosaic (microns)
        eccentricityMicrons;

        % Size of the mRGC mosaic [width, height]
        sizeDegs;
        
        % Either 'right eye' or 'left eye'
        whichEye;

        % The input cone mosaic (@cMosaic)
        inputConeMosaic;

        % Custom degs <-> mm conversions
        customDegsToMMsConversionFunction;
        customMMsToDegsConversionFunction;
        
    end

    % Private properties
    properties (GetAccess=private, SetAccess=private)
        % Size (in degs) of source lattice from which to crop positions for
        % the desired eccentricity
        sourceLatticeSizeDegs = 58;

         % Whether to use parfor in computations
        useParfor = true;

        % Min and max cone positions
        minRFpositionMicrons;
        maxRFpositionMicrons;
    end

    % Public methods
    methods
        
        % Constructor
        function obj = mRGCMosaic(varargin)
            % Parse input
            p = inputParser;
            p.addParameter('name', 'mRGC mosaic', @ischar);
            p.addParameter('centerPolarity', 'ON', @(x)(ischar(x) && (ismember(x, {'ON', 'OFF'}))));
            p.addParameter('whichEye', 'right eye', @(x)(ischar(x) && (ismember(x, {'left eye', 'right eye'}))));
            p.addParameter('eccentricityDegs', [0 0], @(x)(isnumeric(x) && (numel(x) == 2)));
            p.addParameter('sizeDegs', [0.4 0.4], @(x)(isnumeric(x) && (numel(x) == 2)));
            p.addParameter('computeMeshFromScratch', false, @islogical);
            p.addParameter('maxMeshIterations', 100, @(x)(isempty(x) || isscalar(x)));
            p.addParameter('visualizeMeshConvergence', false, @islogical);
            p.addParameter('exportMeshConvergenceHistoryToFile', false, @islogical);
            p.addParameter('inputConeMosaic', [], @(x)(isempty(x) || isa(x, 'cMosaic')));
            p.addParameter('useParfor', true, @islogical);
            p.addParameter('randomSeed', [], @isscalar);
            p.parse(varargin{:});

            obj.name = p.Results.name;
            obj.eccentricityDegs = p.Results.eccentricityDegs;
            obj.sizeDegs = p.Results.sizeDegs;
            obj.randomSeed = p.Results.randomSeed;

            % Parallel computations
            obj.useParfor = p.Results.useParfor;
            

            % Set the inputConeMosaic
            if (isempty(p.Results.inputConeMosaic))
                % An inputConeMosaic was not passed during instantiation,
                % so generate one here. The size of the generated cone
                % mosaic is a little bit larger than the size of the mRGC
                % mosaic so as to allow even the most peripheral mRGCs
                % to have inputs from cones in their surround regions
                maxEccentricityDegs = sqrt(sum((obj.eccentricityDegs + 0.5*obj.sizeDegs).^2,2));
                maxSurroundSigmaDegs = 1/sqrt(2.0) * RGCmodels.CronerKaplan.constants.surroundCharacteristicRadiusFromFitToPandMcells(maxEccentricityDegs);
                % At 2.145 x sigma, a 2D Gaussian is at 1% of the peak
                extraDegsForRGCSurround = 2.146 * maxSurroundSigmaDegs;

                % Compute the input cone mosaic
                obj.inputConeMosaic = cMosaic(...
                    'whichEye', p.Results.whichEye, ...
                    'sizeDegs', obj.sizeDegs + 2 * extraDegsForRGCSurround, ...
                    'eccentricityDegs', obj.eccentricityDegs, ...
                    'opticalImagePositionDegs', 'mosaic-centered');

                  % Watson's non-linear transformations
                  obj.customDegsToMMsConversionFunction = @(x) RGCmodels.Watson.convert.rhoDegsToMMs(x);
                  obj.customMMsToDegsConversionFunction = @(x) RGCmodels.Watson.convert.rhoMMsToDegs(x);

            else
                % An input cone mosaic was passed during instantiation
                % Set the input cone mosaic property to it
                obj.inputConeMosaic = p.Results.inputConeMosaic;
                
                % Match the eccentricity of the mRGCMosaic to that of the input cone mosaic
                oldEccentricityDegs = obj.eccentricityDegs;
                obj.eccentricityDegs = obj.inputConeMosaic.eccentricityDegs;
                fprintf(2, 'mRGCmosaic eccentricity changed from (%2.1f, %2.1f) degs to (%2.1f, %2.1f) degs to interface with input cone mosaic\n', ...
                    oldEccentricityDegs(1), oldEccentricityDegs(2), obj.eccentricityDegs(1), obj.eccentricityDegs(2));

                % Match mm <-> degs conversion functions to those of the inputConeMosaic
                if (~isempty(obj.inputConeMosaic.micronsPerDegreeApproximation))
                    % Generate conversion functions based on the inputConeMosaic.micronsPerDegreeApproximation
                    obj.customDegsToMMsConversionFunction = @(x)(x*1e-3 * obj.inputConeMosaic.micronsPerDegreeApproximation);
                    obj.customMMsToDegsConversionFunction = @(x)(x*1e+3 / obj.inputConeMosaic.micronsPerDegreeApproximation);
                else
                    if (~isempty(obj.inputConeMosaic.customDegsToMMsConversionFunction))
                        % Use the cMosaic conversion functions
                        obj.customDegsToMMsConversionFunction = obj.inputConeMosaic.customDegsToMMsConversionFunction;
                        obj.customMMsToDegsConversionFunction = obj.inputConeMosaic.customMMsToDegsConversionFunction;
                    else
                        % Watson's non-linear transformations
                        obj.customDegsToMMsConversionFunction = @(x) RGCmodels.Watson.convert.rhoDegsToMMs(x);
                        obj.customMMsToDegsConversionFunction = @(x) RGCmodels.Watson.convert.rhoMMsToDegs(x);
                    end
                end

                % Since we were passed an input cone mosaic, we must ensure
                % that the size of the mRGC mosaic is a little bit smaller
                % than the cone mosaic, so as to allow even the most peripheral mRGCs
                % to have inputs from cones in their surround regions
                maxEccentricityDegs = sqrt(sum((obj.eccentricityDegs + 0.5*obj.inputConeMosaic.sizeDegs).^2,2));
                maxSurroundSigmaDegs = 1/sqrt(2.0) * RGCmodels.CronerKaplan.constants.surroundCharacteristicRadiusFromFitToPandMcells(maxEccentricityDegs);
                % At 2.145 x sigma, a 2D Gaussian is at 1% of the peak
                extraDegsForRGCSurround = 2.146 * maxSurroundSigmaDegs;

                oldSizeDegs = obj.sizeDegs;
                obj.sizeDegs = obj.inputConeMosaic.sizeDegs - 2 * extraDegsForRGCSurround;
                fprintf(2, 'mRGCmosaic size changed from %2.1f degs to %2.1f degs to interface with input cone mosaic\n', ...
                    max(oldSizeDegs), max(obj.sizeDegs));
            end

            % Set the whichEye from the input cone mosaic
            obj.whichEye = obj.inputConeMosaic.whichEye;

            % Compute mRGC RF positions
            if (p.Results.computeMeshFromScratch)
                % Re-generate lattice
                obj.regenerateRFPositions(...
                        p.Results.maxMeshIterations,  ...
                        p.Results.visualizeMeshConvergence, ...
                        p.Results.exportMeshConvergenceHistoryToFile);
            else
                % Import positions by cropping a large pre-computed patch
                obj.initializeRFPositions();
            end

            % Remove RFs within the optic disk
            obj.removeRFsWithinOpticNerveHead();

            % Set random seed
            if (isempty(obj.randomSeed))
               rng('shuffle');
            else
               rng(obj.randomSeed);
            end

            % Compute ecc in microns
            obj.eccentricityMicrons = 1e3 * obj.customDegsToMMsConversionFunction(obj.eccentricityDegs);

            % Compute min and max cone position
            obj.minRFpositionMicrons = squeeze(min(obj.rgcRFpositionsMicrons,[],1));
            obj.maxRFpositionMicrons = squeeze(max(obj.rgcRFpositionsMicrons,[],1));
            
            % Wire cones to the RGC RF centers
            obj.wireConesToRFcenters();
        end % Constructor

        % Convert a retinal size in visual degrees at given ecc (also in degs) to its
        % equivalent linear size in retinal microns
        sizeMicrons = sizeVisualDegsToSizeRetinalMicrons(obj,sizeVisualDegs, eccDegs);

        % Convert a retinal size in microns at given ecc (also in microns) to its
        % equivalent angular size in visual degrees
        sizeDegs = sizeRetinalMicronsToSizeVisualDegs(obj, sizeRetinalMicrons, eccMicrons);
    
        % Return the indices of mRGC RFs within a geometry struct appropriate for @regionOfInterest
        rfIndices = indicesOfRFsWithinROI(obj, geometryStruct);

    end % Public methods

    methods (Access=private)
        % Initialize RF positions by importing them from a large previously-computed mesh
        initializeRFPositions(obj);
        
        % Initialize RF positions by regenerating a new mesh. Can be slow.
        regenerateRFPositions(obj, maxIterations, visualizeConvergence, exportHistoryToFile);
        
        % Remove RFs located within the optic disk
        removeRFsWithinOpticNerveHead(obj);

        % Wire cones to the RGC RF centers
        wireConesToRFcenters(obj);

    end % Private methods

end
