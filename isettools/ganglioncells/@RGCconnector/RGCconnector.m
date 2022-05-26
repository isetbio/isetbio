classdef RGCconnector < handle
    % A class to handle connection of a lattice of RGCs to a cMosaic
    %
    % Syntax:
    %   theInputConeMosaic = cMosaic();
    %   rc = RGCconnector(theInputConeMosaic);
    %
    %   theConeMosaic = cMosaic('pigment', pp);
    %   theConeMosaic = cMosaic('sizeDegs',[5,5],'eccentricityDegs',[0,0]);
    %
    % Description
    %

     % Public properties
    properties  (GetAccess=public, SetAccess=public)
        
    end % Public properties

    % Read-only properties
    properties (GetAccess=public, SetAccess=private)
        % The input cone mosaic
        inputConeMosaic;

        % RGC RF positions in microns
        RGCRFpositionsMicrons;

        % RGC RF spacings in microns
        RGCRFspacingsMicrons;

        % Compute struct for computing local cone-to-RGC density ratios
        coneToRGCDensityRatioComputeStruct;

        % Sparse [conesNum x rgcsNum] sparse  connectivity matrix 
        % of cone -> rgc. 
        % To find which cones are connected to a target RGC:
        %  connectivityVector = full(squeeze(obj.coneConnectivityMatrix(:, targetRGC)));
        %  inputConeIDs = find(connectivityVector > 0.01);
        coneConnectivityMatrix = [];

        % Centroids of RGC RFs based on the current cone inputs
        RGCRFcentroidsFromInputs;

        % Params defining wiring pareferences
        wiringParams;

    end

    properties (Constant)
        validRGClatticeModels = {...
            'Watson-midgetRGC' ...
            'Watson-parasolRGC' ...
            };

        defaultWiringParams = struct(...
                'chromaticSpatialVarianceTradeoff', 1.0, ...     % [0: minimize chromatic variance 1: minimize spatial variance]
                'spatialVarianceMetric', 'spatial variance', ... % choose between {'maximal interinput distance', 'spatial variance'}
                'maxNeighborsNum', 6, ...
                'maxNumberOfConesToSwap', 3, ...
                'maxSwapPassesNum', 10, ...
                'maxNeighborNormDistance', 1.5 ...                % max distance to search for neighbors
        );

    end % Constant properties


    % Public methods
    methods
        
        % Constructor
        function obj = RGCconnector(inputConeMosaic, varargin)
            % Validate inputConeMosaic
            assert(isa(inputConeMosaic, 'cMosaic'), 'inputConeMosaic must be a @cMosaic object');
            obj.inputConeMosaic = inputConeMosaic;
            
            % Parse optional input
            p = inputParser;
            p.addParameter('RGCRFpositionsMicrons', [], @(x)((isempty(x)) || (isnumeric(x)&&(size(x,2)==2))));
            p.addParameter('coneToRGCDensityRatio', [], @(x)((isempty(x)) || isnumeric(x)));
            p.addParameter('chromaticSpatialVarianceTradeoff', RGCconnector.defaultWiringParams.chromaticSpatialVarianceTradeoff, @(x)(isscalar(x)&&(x>=0)&&(x<=1)));
            p.addParameter('maxNeighborNormDistance', RGCconnector.defaultWiringParams.maxNeighborNormDistance, @isscalar);
            p.addParameter('maxNumberOfConesToSwap', RGCconnector.defaultWiringParams.maxNumberOfConesToSwap,@(x)(isscalar(x)&&(x>=1)));
            p.addParameter('maxSwapPassesNum', RGCconnector.defaultWiringParams.maxSwapPassesNum, @(x)(isscalar(x)&&(x>=1)));
            p.addParameter('visualizeIntermediateConnectivityStages', false, @islogical);
            p.parse(varargin{:});
            
            RGCRFposMicrons = p.Results.RGCRFpositionsMicrons;
            coneToRGCDensityRatio = p.Results.coneToRGCDensityRatio;
            visualizeIntermediateConnectivityStages = p.Results.visualizeIntermediateConnectivityStages;

            % Update wiringParams struct
            obj.wiringParams = RGCconnector.defaultWiringParams;
            obj.wiringParams.chromaticSpatialVarianceTradeoff = p.Results.chromaticSpatialVarianceTradeoff;
            obj.wiringParams.maxNeighborNormDistance = p.Results.maxNeighborNormDistance;
            obj.wiringParams.maxSwapPassesNum = p.Results.maxSwapPassesNum;
            obj.wiringParams.maxNumberOfConesToSwap = p.Results.maxNumberOfConesToSwap;

            if (isempty(RGCRFposMicrons)) && (isempty(coneToRGCDensityRatio))
                % Nothing was specified, so we initialize with a precomputed RGC lattice (Watson's model)
                modelRGC = 'Watson-midgetRGC';
                fprintf('Instantiating using ''%s'' mRGC lattice\n', modelRGC);
                RGCRFposMicrons = obj.initializeWithPrecomputedLattice(modelRGC);

            elseif (isempty(RGCRFposMicrons)) && (~isempty(coneToRGCDensityRatio))
                % Only a density was specified, so we initialize with a perfect hex RGCRF lattice
                RGCRFposMicrons = obj.initializeWithPerfectHexLattice(coneToRGCDensityRatio);

            elseif (~isempty(RGCRFposMicrons)) && (isempty(coneToRGCDensityRatio))
                % Perhaps issue a warning if there is no overlap with the cMosaic
                fprintf('Instantiating using supplied [%d x 2] lattice of RGC positions\n', ...
                    size(RGCRFposMicrons,1));
            else
                % User specified both a density and a lattice of
                % RGCRFpositions, so raise an error 
                error('@RGCconnector cannot be instantiated both with a set of RGCRFpositions AND a coneToRGCdensityRatio.');
            end

            % Before cropping compute the local cone to RGC density struct
            samplingIntervalMicrons = 3;
            obj.computeConeToRGCDensityRatioComputeStruct(RGCRFposMicrons, samplingIntervalMicrons);

            % Crop positions to lie within the inputConeMosaic
            obj.cropLattice(RGCRFposMicrons);
            
            % Visualize the input mosaics
            if (1==2)
                obj.visualizeInputMosaics();

                % Visualize the cone-to-RGC code density map
                obj.visualizeInputConeToRGCDensityMap();
            end

            % Visualize effective lattice and cone to RGC density map
            if (visualizeIntermediateConnectivityStages)
                obj.visualizeEffectiveConeToRGCDensityMap(900);
            end
            


            % STEP1. Connect cones based on local density.
            obj.connectRGCsToConesBasedOnLocalDensities();

            % Visualize current connectivity
            if (visualizeIntermediateConnectivityStages)
                obj.visualizeCurrentConnectivityState(1001);
            end



            % STEP2. Connect unconnected cones to nearby RGCs
            obj.connectUnconnectedConesToNearbyRGCs(...
                'generateProgressVideo', false);

            % Visualize current connectivity
            if (visualizeIntermediateConnectivityStages)
                obj.visualizeCurrentConnectivityState(1002);
            end


            % STEP3. Transfer 1 cone from RGC1 to a neirboring RGC (RGC2)
            % where the RGC1 has at least N+2inputs, where N = # of inputs in 
            % RGC2, so as to minimize the combined cost.
            % This is the first stage where we fine tune the wiring using
            % params set in the user-supplied wiringParams struct
            obj.transferConesBetweenNearbyRGCsWithUnbalancedInputNumerosities(...
                'generateProgressVideo', ~true);

            % Visualize current connectivity
            if (visualizeIntermediateConnectivityStages)
                obj.visualizeCurrentConnectivityState(1003);
            end


            % STEP 4. Transfer half of the cones from the most populous multi-input RGCs
            % to zero-input RGCs. Here we are utilizing RGCs with 0 inputs
            obj.transferConesFromMultiInputRGCsToZeroInputRGCs(...
                'optimizationCenter', 'visualFieldCenter', ...
                'generateProgressVideo', ~true);

            % Visualize current connectivity
            if (visualizeIntermediateConnectivityStages)
                obj.visualizeCurrentConnectivityState(1004);
            end


            % STEP 5. Swap 1 or more cones  of RGC1 with the same # of cones 
            % in a neighboring RGC2 so at to minimize the combined cost
            obj.swapConesBetweenNearbyRGCs(...
                'optimizationCenter', 'visualFieldCenter', ...
                'generateProgressVideo', ~true);

            % Visualize current connectivity
            if (visualizeIntermediateConnectivityStages)
                obj.visualizeCurrentConnectivityState(1005);
            end


            % Final step of non-overlapping wiring. Remove RGCs on the edges of the patch
            obj.removeRGCsOnPatchPerimeter();

            % Visualize current connectivity
            if (visualizeIntermediateConnectivityStages)
                obj.visualizeCurrentConnectivityState(1009);
            end

        end % Constructor

        % Visualization of input cone mosaics (before any connections are made)
        [hFig, ax, XLims, YLims] = visualizeInputMosaics(obj, varargin);

        % Visualize the current connections between the 2 mosaics
        hFig = visualizeCurrentConnectivityState(obj, figNo);
                
    end % Public methods


    methods (Access=private)
        % Initialize RGCRF positions by importing them from a large
        % previously-computed mesh based on some RGC model
        RGCRFposMicrons = initializeWithPrecomputedLattice(obj, modelRGC);

        % Initialize RGCRF positions with a perfect hexagonal lattice
        % with a desired coneToRGCDensityRatio
        RGCRFposMicrons = initializeWithPerfectHexLattice(obj, coneToRGCDensityRatio);

        % Crop the imported RGCRFlattice to lie within the inputConeMosaic
        cropLattice(obj, RGCRFposMicrons);

        % Compute the obj.coneToRGCDensityRatioComputeStruct. This is done once, just before cropping.
        computeConeToRGCDensityRatioComputeStruct(obj, RGCRFposMicrons, samplingIntervalMicrons);

        % Employ the obj.coneToRGCDensityRatioComputeStruct to compute the
        % local cone-to-RGC density ratios at the current RGC RF positions
        densityRatiosMap = coneToRGCDensityRatiosMap(obj);

        % Find indices of nearby RGCs
        nearbyRGCindices = neihboringRGCindices(obj, theRGCindex);

        % STEP1. Connect RGCs to cones strictly based on local cone-RGC densities
        connectRGCsToConesBasedOnLocalDensities(obj);

        % STEP2. Connect unconnected cones to nearby RGCs
        connectUnconnectedConesToNearbyRGCs(obj, varargin);

        % STEP3. Reassign cone input in nearby RGCs with unbalanced inputs
        transferConesBetweenNearbyRGCsWithUnbalancedInputNumerosities(obj, varargin);

        % STEP4. Transfer half of the cones from the most populous multi-input RGCs to zero-input RGCs
        transferConesFromMultiInputRGCsToZeroInputRGCs(obj, varargin);

        % STEP 5 Swap 1 or more cones  of RGC1 with the same # of cones 
        % in a neighboring RGC2 so at to minimize the combined cost
        swapConesBetweenNearbyRGCs(obj, varargin);

        % Final step of non-overlapping wiring. Remove RGCs on the edges of the patch
        removeRGCsOnPatchPerimeter(obj);

        % Optimize how many and which of theSourceRGCinputConeIndices will
        % be transfered to one of the neighboringRGCindices
        optimizeTransferOfConeInputs(obj, ...
           theSourceRGCindex, theSourceRGCinputConeIndices, theSourceRGCinputConeWeights, ...
           theNeighboringRGCindices, allNeighboringRGCsInputConeIndices, allNeighboringRGCsInputConeWeights);

        % Optimize how many and which of theSourceRGCinputConeIndices will
        % be transfered to a zero input RGC
        optimizeTransferOfConeInputsToZeroInputRGC(obj,...
             theSourceRGCindex, theSourceRGCinputConeIndices, theSourceRGCinputConeWeights, destinationRGCindex);
        
        % Optimize how many and which of theSourceRGCinputConeIndices will
        % be swapped with cones to one of the neighboringRGCindices
        beneficialSwapWasFound = optimizeSwappingOfConeInputs(obj, ...
            theSourceRGCindex, theSourceRGCinputConeIndices, theSourceRGCinputConeWeights, ...
            theNeighboringRGCindices, allNeighboringRGCsInputConeIndices, allNeighboringRGCsInputConeWeights);

        % Update the connectivityMatrix, by disconnecting
        %   indexOfConeToBeReassigned  FROM  sourceRGCindex
        % and connecting 
        %   indexOfConeToBeReassigned  to the destinationRGCindex
        transferConeFromSourceRGCToDestinationRGC(obj, ...
             indexOfConeToBeReassigned, sourceRGCIndex, destinationRGCindex);

        % Update the connectivityMatrix by swapping cones
        %       sourceRGCconeIndicesToBeSwapped from sourceRGCindex
        % with cones
        %       destinationRGCconeIndicesToBeSwapped form destinationRGCindex
        swapConesFromSourceRGCWithConesOfDestinationRGC(obj, ...
            sourceRGCconeIndicesToBeSwapped, sourceRGCindex, ...
            destinationRGCconeIndicesToBeSwapped, destinationRGCindex)
    

        % Compute the cost for an RGC to maintain its cone inputs
        [cost, spatialCostComponent, chromaticCostComponent] = ...
            costToMaintainInputs(obj, inputConeIndices, inputConeWeights, targetRGCSpacingMicrons);

        % Update the centroids of all RGCs in the RGClist
        updateCentroidsFromInputs(obj, RGClist);

        % Update RGCRFspacings for all RGCs beased on their current centroids
        updateRGCRFspacingsBasedOnCurrentCentroids(obj);

        % Visualize the cones of the input cone mosaic using a custom shape
        % cone outline
        visualizeConePositions(obj, ax, shapeOutline, varargin);

        % Visualize the input cones to each RGC
        visualizeRGCinputs(obj, ax, varargin);

        % Visualization of the cone-to-RGC density ratios at the positions of the RGCs
        visualizeEffectiveConeToRGCDensityMap(obj, figNo);

        % Visualization of the continuous cone-to-RGC density map
        visualizeInputConeToRGCDensityMap(obj);

        % Visualization of the connectivity between cones and RGCRFs
        [hFig, ax, XLims, YLims] = visualizeConnectivity(obj, varargin);

        % Visualize spatial variance cost statistics
        visualizeSpatialVarianceCostStatistics(obj, axSpatial, spatialVarianceCost);

        % Visualize chromatic variance cost statistics
        visualizeChromaticVarianceCostStatistics(obj, axChromatic, chromaticVarianceCost);
    end


    methods (Static)
        % Compute methods
        [D,idx] = pdist2(A, B, varargin);
        d = maximalInterInputDistance(coneRFpos);
    
        % Indices of points are inside & on the boundary defined by a select subset of points
        [insideBoundaryPointIndices, onBoundaryPointIndices] = ...
            pointsInsideBoundaryDefinedBySelectedPoints(allPointPositions, selectedPointIndices);

        % Visualization methods
        transparentContourPlot(axesHandle, spatialSupportXY, zData, ...
          zLevels, faceAlpha, cmap, lineStyle, lineWidth);

        [f,v] = facesAndVertices(positions, spacings, diskOutline);
    end
end

