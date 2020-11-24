classdef mRGCmosaic < handle
% Create a midget RGC mosaic connected to a cone mosaic

    properties (Constant)
        LCONE_ID = 2;
        MCONE_ID = 3;
        SCONE_ID = 4;
    end
    
    % Public properties
    properties
        % 'noiseFlag'     - String. Default 'random'. Add Gaussian noise
        %                 (default) or not. Valid values are {'random', 'frozen', 'none'}.
        noiseFlag = 'random';
        
        % 'noiseFactor'   - Float. Default 0.2. Sigma of Gaussian noise = noiseFactor * max(noiseFreeResponse)
        noiseFactor = 0.2;
    end
    
    % Read-only properties
    properties (GetAccess=public, SetAccess=private)
        
        % The eccentricity of the mosaic, in degrees
        eccentricityDegs;
        
        % The size of the mosaic, in degrees
        sizeDegs;
        
        % Eye, left or right
        whichEye;
        
        % The input cone mosaic
        inputConeMosaic;
        
        % Sparse matrix [nCones x mRGC] storing the exclusive connections
        % between the n-th cone to m-th RGC center subregion 
        % (1==connected, 0==disconencted)
        coneConnectivityMatrix;
        
        % Struct containing sparse matrices with weights of cone connections
        % to the RGC center & surround subregions
        coneWeights;
        
        % [m x 2] matrix of RGC positions, in microns
        rgcRFpositionsMicrons;
        
        % [m x 2] matrix of RGC positions, in degrees
        rgcRFpositionsDegs;
    end
    
    % Private properties
    properties (GetAccess=private, SetAccess=private)
        % The metadata of the input cone mosaic (used for connecting cones 
        % to the subregions of the mRGC cell receptive fields.
        inputConeMosaicMetaData;
        
        % Size (in degs) of source lattice from which to crop positions for
        % the desired eccentricity
        sourceLatticeSizeDegs = 45;
        
        % [m x 1] matrix of local spacing for each RGC, in microns
        rgcRFspacingsMicrons;
        
        % [m x 1] matrix of local spacing for each RGC, in degs
        rgcRFspacingsDegs;
        
        % Synthesized RGC RF params
        synthesizedRFparams;
        
        % Imported cone and mRGC positions
        importedData;
    end
    
    % Public methods
    methods
        % Constructor
        function obj = mRGCmosaic(eccentricityDegs, sizeDegs, whichEye, varargin)
            
            p = inputParser;
            p.addRequired('eccentricityDegs', @(x)(isnumeric(x) && (numel(x) == 2)));
            p.addRequired('sizeDegs', @(x)(isnumeric(x) && (numel(x) == 2)));
            p.addRequired('whichEye', @(x)(ischar(x) && (ismember(x, {'left', 'right'}))));
            p.addParameter('viewTesselationMaps', false, @islogical);
            p.addParameter('coneSpecificityLevel', 100, @(x)(isscalar(x) && (x>=0) && (x<=100)));
            
            % Allow unmatched key-value pairs. These are to be used by
            % cone mosaic.
            p.KeepUnmatched = true;
            p.parse(eccentricityDegs, sizeDegs, whichEye, varargin{:});
    
            % Make sure any unmatched key-value pairs are valid
            unmatchedParameterNames = fieldnames(p.Unmatched);
            validUnmatchedParameterNames = {...
                'coneMosaicResamplingFactor' ...
                'coneMosaicSpatialDensity' ...
                'coneMosaicIntegrationTime' };
            
            for iUnmatched = 1:numel(unmatchedParameterNames)
                assert(ismember(unmatchedParameterNames{iUnmatched}, validUnmatchedParameterNames), ...
                    fprintf('Parameter ''%'' is unexpected.', unmatchedParameterNames{iUnmatched}));
            end
            
            % Set properties
            viewTesselationMaps = p.Results.viewTesselationMaps;
            coneSpecificityLevel = p.Results.coneSpecificityLevel;
            
            obj.eccentricityDegs = eccentricityDegs;
            obj.sizeDegs = sizeDegs;
            obj.whichEye = whichEye;
            
            fprintf('Generating input cone mosaic. Please wait ...\n');
                
            % Import cone and mRGC RF positions by cropping large
            % eccentricity-varying cone and mRGC lattices
            [obj.importedData.coneRFpositionsMicrons, ...
             obj.importedData.coneRFpositionsDegs, ...
             obj.importedData.rgcRFpositionsMicrons, ...
             obj.importedData.rgcRFpositionsDegs, ...
             extraDegsForRGCSurround] =  mRGCmosaic.importConeAndRGCpositions(...
                        obj.sourceLatticeSizeDegs, ...
                        eccentricityDegs, ...
                        sizeDegs, ...
                        whichEye);

            % Generate a regular hex mosaic to serve as the
            % input cone mosaic with a custom mean cone spacing 
            % (equal to the median spacing within the imported
            % coneRFpositionsMicrons) and custom optical density and
            % macular pigment appropriate for the eccentricityDegs
            generationMode = 'equivalent regular hex';
            [obj.inputConeMosaic, ...
             obj.inputConeMosaicMetaData] = mRGCmosaic.generateInputConeMosaic(...
                        generationMode, ...
                        eccentricityDegs, ...
                        sizeDegs, ...
                        extraDegsForRGCSurround, ...
                        obj.importedData.coneRFpositionsMicrons, ...
                        varargin{:});
            
            % Wire cones to RGC centers, computing the cone-to-RGC-center
            % connectivity matrix, and the resulting RGC RF positions and
            % spacings
            [obj.coneConnectivityMatrix, ...
             obj.rgcRFpositionsDegs, obj.rgcRFpositionsMicrons, ...
             obj.rgcRFspacingsDegs, obj.rgcRFspacingsMicrons] = mRGCmosaic.wireInputConeMosaicToRGCcenters(...
                obj.importedData.rgcRFpositionsDegs, obj.importedData.rgcRFpositionsMicrons,  ...
                obj.inputConeMosaicMetaData.conePositionsDegs, ...
                obj.inputConeMosaicMetaData.conePositionsMicrons, ...
                obj.inputConeMosaicMetaData.coneSpacingsMicrons, ...
                obj.inputConeMosaicMetaData.coneTypes, ...
                obj.inputConeMosaicMetaData.indicesOfConesNotConnectingToRGCcenters, ...
                coneSpecificityLevel, viewTesselationMaps);
            
            % Compute weights of connections between cones and RGC
            % center/surround subregions, and update tge RGC RF positions
            % based on the connectivity
            [obj.coneWeights, obj.rgcRFpositionsDegs, obj.synthesizedRFparams] = mRGCmosaic.computeConeWeights(...
                obj.inputConeMosaicMetaData.conePositionsDegs, ...
                obj.inputConeMosaicMetaData.coneTypes, ...
                obj.coneConnectivityMatrix);
            
            % Update RGCRF positions and spacings
            obj.rgcRFpositionsMicrons = 1e3*RGCmodels.Watson.convert.rhoDegsToMMs(obj.rgcRFpositionsDegs);
            obj.rgcRFspacingsDegs = RGCmodels.Watson.convert.positionsToSpacings(obj.rgcRFpositionsDegs);
            obj.rgcRFspacingsMicrons = RGCmodels.Watson.convert.positionsToSpacings(obj.rgcRFpositionsMicrons);
        end
        
        % Method to compute the response of the RGC mosaic for the passed cone mosaic response
        [mRGCresponses, temporalSupportSeconds] = compute(obj, coneMosaicResponses, timeAxis, varargin);
        
        % Visualize mRGC positions together with imported ecc-varying cone positions
        % and together with cone positions in the equivalent employed reg-hex mosaic
        visualizeInputPositions(obj);
        
        % Method to visualize the tesselation of the input cone mosaic by
        % the RF centers of the RGC mosaic
        visualizeConeMosaicTesselation(obj, domain);
        
        % Method to visualize the synthesized RF params, retinal and visual
        visualizeSynthesizedParams(obj);
        
        % Method to visualize the cone weights to each RGC
        visualizeConeWeights(obj);
        
        % Method to visualize the mosaic activation
        visualizeActivationMap(obj, axesHandle, activation, varargin);
        
        % Method to visualize mRGC mosaic responses
        visualizeResponses(obj, temporalSupportSeconds, mRGCMosaicResponse, varargin);
        
        % Method to visualize mRGC responses in a matrix matching the RGC positions
        visualizeResponseMatrix(obj, independentVariable, responseMatrix, varargin);
        
        visualizeCorrespondenceBetweenMappedAndSynthesizedModelParams(obj, visuallyMappedModelParams);
        
        % Method to generate axes positions for visualizing mRGC responses in a matrix matching the RGC positions 
        axPos = axesMatrixPosition(obj);
    end
    
    % Static methods
    methods (Static)
        % Static method to import cone and mRGC positions from pre-computed
        % lattices that have the desired size and are centered at the desired eccentricity in the desired eye
        [coneRFpositionsMicrons, coneRFpositionsDegs, ...
         rgcRFpositionsMicrons,  rgcRFpositionsDegs, extraDegsForRGCSurround] = ...
            importConeAndRGCpositions(sourceLatticeSizeDegs, eccentricityDegs, sizeDegs, whichEye);
        
        % Static method to generate a cone mosaic from the imported cone positions
        [theConeMosaic, theConeMosaicMetaData] = ...
            generateInputConeMosaic(generationMode, eccentricityDegs, sizeDegs, ...
            extraDegsForRGCSurround, coneRFpositionsMicrons, varargin); 
        
        % Static method to wire cones to the RGC RF centers
        [connectivityMatrix, rgcRFpositionsDegs, rgcRFpositionsMicrons, rgcRFspacingsDegs, rgcRFspacingsMicrons] = ...
            wireInputConeMosaicToRGCcenters(rgcRFpositionsDegs, rgcRFpositionsMicrons, ...
            conePositionsDegs, conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
            indicesOfConesNotConnectingToRGCcenters, coneSpecificityLevel, viewTesselationMaps);
        
        % Static method to compute weights of connections b/n cones and
        % center-surround RF subregions. Also update RGC RF positions
        [coneWeights, rgcPositionsDegsFromConnectivity, synthesizedRFparams] = computeConeWeights(conePositionsDegs, ...
            coneTypes, connectivityMatrix, eccentricityDegs, sizeDegs);
    
        % Static method to render the map of how RGC centers tesellate cones
        renderTesselationMap(figNo, axesHandle, coneRFpositions, coneRFspacings, coneTypes, rgcRFpositions, rgcRFspacings, coneConnectivityMatrix, domain);
    end
end

