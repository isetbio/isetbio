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

        % Locally-averaged RGCRFspacingsMicrons;
        localRGCRFspacingsMicrons;

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

    % Constant properties
    properties (Constant)
        validRGClatticeModels = {...
            'Watson-midgetRGC' ...
            'Watson-parasolRGC' ...
            };

        defaultWiringParams = struct(...
                'chromaticSpatialVarianceTradeoff', 1.0, ...          % [0: minimize chromatic variance 1: minimize spatial variance]
                'rfCentroidOverlapPenaltyFactor', 1, ...                      % Penalty for overlapping centroids
                'RcToRGCseparationRatio', 1.0, ...                    % overlap of RFs (1 = no overlap)
                'spatialVarianceMetric', 'spatial variance', ...      % choose between {'maximal interinput distance', 'spatial variance'}
                'maxNeighborsNum', 6, ...                             % max numboer of neighboring RGCs
                'maxNeighborNormDistance', 1.5, ...                   % max distance to search for neighboring RGCs
                'maxMeanConeInputsPerRGCToConsiderSwapping', 10, ...  % Do cone swapping only if the mean cones/RGC less than or equal to this number
                'maxNumberOfConesToSwap', 6, ...                      % Only swap up to this many cones
                'maxPassesNum', 50 ...     
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
            p.addParameter('RGCRFspacingsMicrons', [], @(x)((isempty(x)) || (isnumeric(x))));
            p.addParameter('coneToRGCDensityRatio', [], @(x)((isempty(x)) || isnumeric(x)));
            p.addParameter('RcToRGCseparationRatio', 0, @(x)(isscalar(x)&&(x>=1)));
            p.addParameter('chromaticSpatialVarianceTradeoff', RGCconnector.defaultWiringParams.chromaticSpatialVarianceTradeoff, @(x)(isscalar(x)&&(x>=0)&&(x<=1)));
            p.addParameter('rfCentroidOverlapPenaltyFactor', 1, @(x)(isscalar(x)&&(x>=0)));
            p.addParameter('maxNeighborNormDistance', RGCconnector.defaultWiringParams.maxNeighborNormDistance, @isscalar);
            p.addParameter('maxNumberOfConesToSwap', RGCconnector.defaultWiringParams.maxNumberOfConesToSwap,@(x)(isscalar(x)&&(x>=1)));
            p.addParameter('maxMeanConeInputsPerRGCToConsiderSwapping', RGCconnector.defaultWiringParams.maxMeanConeInputsPerRGCToConsiderSwapping, @(x)(isscalar(x)&&(x>=1)));
            p.addParameter('maxPassesNum', RGCconnector.defaultWiringParams.maxPassesNum, @(x)(isscalar(x)&&(x>=1)));
            p.addParameter('visualizeIntermediateConnectivityStages', false, @islogical);
            p.parse(varargin{:});
            
            RGCRFposMicrons = p.Results.RGCRFpositionsMicrons;
            RGCRFspacingsMicrons = p.Results.RGCRFspacingsMicrons;
            coneToRGCDensityRatio = p.Results.coneToRGCDensityRatio;
            
            visualizeIntermediateConnectivityStages = p.Results.visualizeIntermediateConnectivityStages;

            % Update wiringParams struct
            obj.wiringParams = RGCconnector.defaultWiringParams;
            obj.wiringParams.chromaticSpatialVarianceTradeoff = p.Results.chromaticSpatialVarianceTradeoff;
            obj.wiringParams.rfCentroidOverlapPenaltyFactor = p.Results.rfCentroidOverlapPenaltyFactor;
            obj.wiringParams.maxNeighborNormDistance = p.Results.maxNeighborNormDistance;
            obj.wiringParams.maxPassesNum = p.Results.maxPassesNum;
            obj.wiringParams.maxNumberOfConesToSwap = p.Results.maxNumberOfConesToSwap;
            obj.wiringParams.maxMeanConeInputsPerRGCToConsiderSwapping = p.Results.maxMeanConeInputsPerRGCToConsiderSwapping;
            obj.wiringParams.RcToRGCseparationRatio = p.Results.RcToRGCseparationRatio;
            
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
            obj.computeConeToRGCDensityRatioComputeStruct(RGCRFposMicrons, RGCRFspacingsMicrons,samplingIntervalMicrons);

            % Crop positions to lie within the inputConeMosaic
            obj.cropLattice(RGCRFposMicrons);
            
            if (1==2)
                [hFig, ax, XLims, YLims] = obj.visualizeInputMosaics();
                obj.visualizeConnectivity('figureHandle', hFig, 'axesHandle', ax, 'XLims', XLims, 'YLims', YLims);
                set(hFig, 'Position', [ 10 10 1030 700]);
                NicePlot.exportFigToPDF('mosaics.pdf', hFig, 300);
            end
            
            
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
            

            % STEP3. Transfer cones (1 at a time) from RGC1 to a neirboring RGC (RGC2)
            % where the RGC1 has at least N+2inputs, where N = # of inputs in 
            % RGC2, so as to minimize the combined cost for RGC1+RGC2.
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


            % STEP 5. Swap 1 or more cones of RGC1 with the same # of cones 
            % in a neighboring RGC2 so at to minimize the combined cost
            obj.swapConesBetweenNearbyRGCs(...
                'optimizationCenter', 'visualFieldCenter', ...
                'generateProgressVideo', ~true);

            % Visualize current connectivity
            if (visualizeIntermediateConnectivityStages)
                obj.visualizeCurrentConnectivityState(1005);
            end

            % STEP 6. Remove RFs on perimeter
            %obj.removeRGCsOnPatchPerimeter();

            % Allow for overlapping of cone inputs
            obj.divergeConeOutputsToMultipleNearbyRGCs();
            if (visualizeIntermediateConnectivityStages)
                obj.visualizeCurrentConnectivityState(1006);
            end
        end % Constructor

        % Visualization of input cone mosaics (before any connections are made)
        [hFig, ax, XLims, YLims] = visualizeInputMosaics(obj, varargin);

        % Visualize the current connections between the 2 mosaics
        hFig = visualizeCurrentConnectivityState(obj, figNo, varargin);
         
        % Visualize the full connectivity of a single RGC
        hFig = visualizeConePoolingWithinRFcenter(obj, iRGC, varargin);

        % Visualize the cone pooling & overlap within RGCs nearby an RGC
        [hFig,visualizedNeighborsNum] = visualizeConePoolingWithinNeighboringRGCs(obj, iRGC, varargin);
        
        % Compute input maintenance costs across entire RGC mosaic
        [totalCost, spatialCost, chromaticCost] = computeInputMaintenanceCostAcrossEntireMosaic(obj);

        % Diverge cone outputs to multiple nearby RGCs (creating RF overlap)
        divergeConeOutputsToMultipleNearbyRGCs(obj, varargin);

        % Return indices of nearby RGCs
        nearbyRGCindices = neihboringRGCindices(obj, theRGCindex, varargin);
        

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
        computeConeToRGCDensityRatioComputeStruct(obj, RGCRFposMicrons, RGCRFspacingsMicrons, samplingIntervalMicrons);

        % Employ the obj.coneToRGCDensityRatioComputeStruct to compute the
        % local cone-to-RGC density ratios at the current RGC RF positions
        densityRatiosMap = coneToRGCDensityRatiosMap(obj);


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


        % Remove RGCs on the edges of the patch
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
        [beneficialSwapWasFound, costReduction] = optimizeSwappingOfConeInputs(obj, ...
            theSourceRGCindex, theSourceRGCinputConeIndices, theSourceRGCinputConeWeights, ...
            theNeighboringRGCindices, allNeighboringRGCsInputConeIndices, allNeighboringRGCsInputConeWeights, dryRunOnly);

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
            costToMaintainInputs(obj, inputConeIndices, inputConeWeights, localRGCSpacingMicrons);

        projectedCostFromOverlap = costToMaintainOverlappingInputs(obj, ...
           neighboringRGCindex, neighboringRGCconeIndices, neighboringRGCconeWeights, ...
           sourceRGCindex, sourceRGCconeIndices, sourceRGCconeWeights);

        % Update the centroids of all RGCs in the RGClist
        updateCentroidsFromInputs(obj, RGClist);

        % Update RGCRFspacings for all RGCs beased on their current centroids
        updateLocalRGCRFspacingsBasedOnCurrentCentroids(obj);

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

        % Visualize the pooled cone apertures of a single RGC
        [xSupport, ySupport, rfProfile2D, xTicks, yTicks] = ...
            visualizeConeAperturePooling(obj, iRGC, axConeWiring, axConeApertures, ...
            xSupport, ySupport, xTicks, yTicks, visualizedFieldOfViewMicrons, visualizedConesNum, colorString)

        % Visualize the 2D overlap between 2 RGC RFs
        visualizeRFoverlap2D(obj, iRGC, nearbyRGCindex, ...
            rfProfile2DmainRGC, rfProfile2DnearbyRGC, ...
            xSupport, ySupport, xTicks, yTicks, axRFOverlap2D);
        
        % Visualize spatial variance cost statistics
        visualizeSpatialVarianceCostStatistics(obj, axSpatial, spatialVarianceCost);

        % Visualize chromatic variance cost statistics
        visualizeChromaticVarianceCostStatistics(obj, axChromatic, chromaticVarianceCost);
    end % Private methods


    methods (Static)
        % Compute methods
        [D,idx] = pdist2(A, B, varargin);
        d = maximalInterInputDistance(coneRFpos);
        c = weightedMean(data, weights);
        overlapCoeff = overlap(weights1, weights2);
    
        shiftedPositions();
        
        % Indices of points are inside & on the boundary defined by a select subset of points
        [insideBoundaryPointIndices, onBoundaryPointIndices] = ...
            pointsInsideBoundaryDefinedBySelectedPoints(allPointPositions, selectedPointIndices);

        % Determine whether convergece is achieved
        convergenceIsAchieved = convergenceAchieved(netTotalCostSequence);

        % Visualization methods
        visualizeConvergence(currentPass, netTotalCostInitial, netTotalCost, ...
            netSpatialCostInitial, netSpatialCost, ...
            netChromaticCostInitial, netChromaticCost, ...
            netReassignments, maxPassesNum);

        transparentContourPlot(axesHandle, spatialSupportXY, zData, ...
          zLevels, faceAlpha, cmap, lineStyle, lineWidth);

        shadedAreaPlot(ax,x,y, baseline, faceColor, edgeColor, faceAlpha, lineWidth, lineStyle)
    
        [f,v] = facesAndVertices(positions, spacings, diskOutline);
    end % Static methods
end
