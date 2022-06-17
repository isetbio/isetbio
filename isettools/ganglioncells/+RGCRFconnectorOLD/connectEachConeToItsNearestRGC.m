function  [RGCRFinputs, RGCRFweights, availableZeroInputRGCsNum] = ...
                connectEachConeToItsNearestRGC(theInputConeMosaic, RGCRFpos, ...
                        localConeToRGCDensityRatioStruct, ...
                        wiringParams, visualizationParams)
% Connect each of a set of input cones (L- or M-) to a set of RGC RF centers
%
% Syntax:
%   [RGCRFinputs, RGCRFweights, availableZeroInputRGCsNum] = ...
%        RGCRFconnector.connectEachConeToItsNearestRGC(...
%                    theInputConeMosaic, RGCRFpos, ...
%                    localConeToRGCDensityRatioStruct, ...
%                    wiringParams, visualizationParams)
%
% Description:
%   Connect each of a set of input cones (L- or M-) to a set of RGC RF centers
%
% Inputs:
%    theInputConeMosaic                 - the input cone mosaic
%    RGCRFposMicrons                    - [M x 2] matrix of (x,y) positions (in microns) of M target RGC RF centers
%    localConeToRGCDensityRatioStruct   - Struct with info to compute the cone/RGC density ratio at any (x,y) position
%    wiringParams                       - Struct with the following wiring params
%       chromaticSpatialVarianceTradeoff   - Chromatic-SpatialVariance tradefoff
%       sequentialConeTypeWiring           - Logical. Whether to wire cone types sequentially or not
%       maxNearbyRGCsNum                   - Scalar. Max number of nearby RGCs to look for assignment, if a cone cannot
%                                         be assigned to its closest RGC because that RGC already has one input cone
%    visualizationParams                 - Struct with the following visualization params
%       generateProgressionVideo           - Boolean. Whether to generate a video showning the different wiring steps
%       exportStagesAsPDF                  - Boolean. Whether to generate a PDF at the end of each wiring stage
%
% Outputs:
%    RGCRFinputs                - Cell array with indices of the input cones, one cell per each target RGC RF center
%    RGCRFweights               - Cell array with weights of the input cones, one cell per each target RGC RF center   
%    availableZeroInputRGCsNum  - Number of the M target RGCs with zero inputs 
%
% Optional key/value pairs
%   none
%   
% History:
%   5/11/2022       NPC     Wrote it
%

    % Unpack inputs
    allConePositions = theInputConeMosaic.coneRFpositionsMicrons;
    allConeSpacings = theInputConeMosaic.coneRFspacingsMicrons;
    allConeTypes = theInputConeMosaic.coneTypes;

    
    % Stage 1A. Connect each RGC to M single cones: the closest ones. Also, keep a track of 
    % cone indices that were not connected.
    connectedConeIndices = [];
    
    if ((visualizationParams.exportStagesAsPDF) || (visualizationParams.generateProgressionVideo))
        % PDFfilename
        if (wiringParams.sequentialConeTypeWiring)
            PDFfilename = sprintf('WiringSteps_W%1.2f_SequentialConeTypeWiring_Stage1', wiringParams.chromaticSpatialVarianceTradeoff);
        else
            PDFfilename = sprintf('WiringSteps_W%1.2f_NonSequentialConeTypeWiring_Stage1', wiringParams.chromaticSpatialVarianceTradeoff);
        end
    
        % Video setup
        if (visualizationParams.generateProgressionVideo)
            % Video setup
            videoOBJ = VideoWriter(PDFfilename, 'MPEG-4');
            videoOBJ.FrameRate = 10;
            videoOBJ.Quality = 100;
            videoOBJ.open();
        else
            videoOBJ = '';
        end

        hFig = figure(100); clf;
        set(hFig, 'Name', 'Stage 1A', 'Color', [1 1 1]);
        ax = subplot('Position', [0.1 0.1 0.8 0.8]);
    end

    % Store the cone inputs in a temporary array
    tmpRGCRFinputs = cell(1,size(RGCRFpos,1));

    % Compute the local code/RGC density ratios for all RGS
    localConeToRGCDensityRatios = RGCRFconnector.localConeToRGCDensityRatiosAtScatteredPositions(...
                    localConeToRGCDensityRatioStruct,...
                    RGCRFpos);

    if (wiringParams.sequentialConeTypeWiring)
        % Connect the M-cones during the first pass
        targetConeTypes = [cMosaic.MCONE_ID];
        [RGCRFpos, tmpRGCRFinputs, connectedConeIndices] = DoIt(...
           RGCRFpos, tmpRGCRFinputs, ...
           localConeToRGCDensityRatios, ...
           connectedConeIndices, ...
           allConePositions,  ...
           allConeTypes, targetConeTypes);
    
        % Connect the L-cones during the second pass
        targetConeTypes = [cMosaic.LCONE_ID];
        [RGCRFpos, tmpRGCRFinputs, connectedConeIndices] = DoIt(...
           RGCRFpos, tmpRGCRFinputs, ...
           localConeToRGCDensityRatios, ...
           connectedConeIndices, ...
           allConePositions,  ...
           allConeTypes, targetConeTypes);

    else
        % Connect both L and M-cones at the same pass
        targetConeTypes = [cMosaic.LCONE_ID cMosaic.MCONE_ID];
        [RGCRFpos, tmpRGCRFinputs, connectedConeIndices] = DoIt(...
           RGCRFpos, tmpRGCRFinputs, ...
           localConeToRGCDensityRatios, ...
           connectedConeIndices, ...
           allConePositions,  ...
           allConeTypes, targetConeTypes);
    end



    % List of unconnected cone indices
    unconnectedConeIndices = setdiff(1:size(allConePositions,1), connectedConeIndices);
    unconnectedConeIndices = setdiff(unconnectedConeIndices, theInputConeMosaic.sConeIndices);

    if ((visualizationParams.exportStagesAsPDF) || (visualizationParams.generateProgressionVideo))
        stageName = 'Stage-1A';
        [hFig, ax] = setupFigure(stageName);
        
        [~, ~, XLims, YLims] = RGCRFconnector.visualizeWiringOfInputConeMosaicToRFcenters(...
            tmpRGCRFinputs, [], ...
            theInputConeMosaic, ...
            'titleString', stageName, ...
            'pauseAfterEachRGCisRendered', false, ...
            'superimposeConeInputWiring', true, ...
            'figureHandle', hFig, ...
            'videoOBJ', videoOBJ, ...
            'axesHandle', ax);

        if (visualizationParams.exportStagesAsPDF)
            NicePlot.exportFigToPDF(sprintf('%s-Stage1A.pdf', PDFfilename), hFig, 300);
        end
        pause
    end

    
    fprintf('%d out of %d L/M cones could not be connected to their closest RGC because the closest RGC already had the expected L/M cone attached\n', ...
        numel(unconnectedConeIndices), numel((allConeTypes ~= cMosaic.SCONE_ID)));
    fprintf('Trying to connect them to nearby RGCs, minimizing the cost\n');
    

    % Sort unconnected cone indices according to their distance from the
    % center of the mosaic
    d = sum((bsxfun(@minus, allConePositions(unconnectedConeIndices,:), theInputConeMosaic.eccentricityMicrons)).^2,2);
    [~,idx] = sort(d, 'ascend');
    unconnectedConeIndices = unconnectedConeIndices(idx);

    % Stage 1B. Deal with the unconnected cone indices
    for iCone = 1:numel(unconnectedConeIndices)
        theConeIndex = unconnectedConeIndices(iCone);
        theConeRFPosition = allConePositions(theConeIndex,:);

        % Find the maxNearbyRGCsNum neigboring RGCs 
        [~, nearestRGCindices] = pdist2(RGCRFpos, theConeRFPosition, '', ...
            'smallest', wiringParams.maxNearbyRGCsNum);

        % Check to see if any of these RGCs have 0 inputs.
        targetRGCindex = [];
        for iRGC = 1:numel(nearestRGCindices)
            % Retrieve this cell's cone inputs
            theRGCindex = nearestRGCindices(iRGC);
            inputConeIndices = tmpRGCRFinputs{theRGCindex};
            if (numel(inputConeIndices) == 0)
                targetRGCindex = theRGCindex;
            end
        end

        if (~isempty(targetRGCindex))
            % Found an RGC with zero cone inputs. Use it
            tmpRGCRFinputs{targetRGCindex} = theConeIndex;

            % Update the position of this RGCRF
            RGCRFpos(targetRGCindex,:) = theConeRFPosition;
                  
        else
            % All near RGC have at least one cone input.
            % Compute the cost of assigning this cone to each of these RGCs
            projectedCosts = zeros(1,numel(nearestRGCindices));
    
            % Compute the theoretical cone/RGC density ratios at the
            % positions of the RGCs
            localConeToRGCDensityRatios = ...
                RGCRFconnector.localConeToRGCDensityRatiosAtScatteredPositions(...
                    localConeToRGCDensityRatioStruct,...
                    RGCRFpos(nearestRGCindices,:));
            

            for iRGC = 1:numel(nearestRGCindices)
                % Retrieve this cell's cone inputs
                theRGCindex = nearestRGCindices(iRGC);
                inputConeIndices = tmpRGCRFinputs{theRGCindex};
    
                % To test the new cost, add this cone as an input
                inputConeIndices(numel(inputConeIndices)+1) = theConeIndex;
                inputConePositions = allConePositions(inputConeIndices,:);
                inputConeSpacings = allConeSpacings(inputConeIndices);
                inputConeTypes = allConeTypes(inputConeIndices);
                inputConeWeights = ones(1, numel(inputConeIndices));
    
                % Compute new cost, had this cone been assigned to this RGC
                projectedCosts(iRGC) = RGCRFconnector.costToMaintainInputCones(...
                    wiringParams.chromaticSpatialVarianceTradeoff, ...
                    inputConePositions, inputConeSpacings, ...
                    inputConeTypes, inputConeWeights, ...
                    localConeToRGCDensityRatios(iRGC));
            end % iRGC

            % Find the RGCs with the min projected cost
            minCost = min(projectedCosts);
            idx = find(projectedCosts == minCost);
            candidateTargetRGCindices = nearestRGCindices(idx);
            
            if (numel(candidateTargetRGCindices)==1)
                targetRGCindex = candidateTargetRGCindices(1);
            else
                % More than one of the nearest RGCs have the same minimal cost
                % Selection strategy
                strategy = 'smallest maximal inter-input distance';
                %strategy = 'nearest RGC';
    
                switch (strategy)
                    case 'smallest maximal inter-input distance'
                        % More than one RGCs have the same minimum cost
                        % Pick the one that has the smallest maximal inter-input distance
                        maximalInterInputDistances = zeros(1,numel(candidateTargetRGCindices));
                        for kkk = 1:numel(candidateTargetRGCindices)
                            % Find the inputs to this RGC
                            theInputConeIndicesForThisRGC = tmpRGCRFinputs{candidateTargetRGCindices(kkk)};
                            d = pdist2(allConePositions(theInputConeIndicesForThisRGC,:), theConeRFPosition);
                            maximalInterInputDistances(kkk) = max(d(:));
                        end
                        [~, kkk] = min(maximalInterInputDistances);
                        targetRGCindex = candidateTargetRGCindices(kkk);
                        
                    case 'nearest RGC'
                        % More than one RGCs have the same minimum cost
                        % Pick the one whose centroid is closest to the cone 
                        d = pdist2(RGCRFpos(candidateTargetRGCindices,:), theConeRFPosition);
                        [~, kkk] = min(d(:));
                        targetRGCindex = candidateTargetRGCindices(kkk);
    
                    otherwise
                        error('Unknown strategy %''%s''.', strategy);
                end
            end

            % Assign theConeIndex to this RGC
            inputConeIndices = tmpRGCRFinputs{targetRGCindex};
            inputConeIndices(numel(inputConeIndices)+1) = theConeIndex;
            tmpRGCRFinputs{targetRGCindex} = inputConeIndices;
    
            % Update the position of this RGCRF
            RGCRFpos(targetRGCindex,:) = RGCRFconnector.centroidsFromConeInputs({tmpRGCRFinputs{targetRGCindex}}, [], allConePositions);
        end

        if (visualizationParams.generateProgressionVideo)
            stageName = 'Stage-1B';
            [hFig, ax] = setupFigure(stageName);
            [~, ~, XLims, YLims] = RGCRFconnector.visualizeWiringOfInputConeMosaicToRFcenters(...
                tmpRGCRFinputs, [], ...
                theInputConeMosaic, ...
                'titleString', stageName, ...
                'superimposeConeInputWiring', true, ...
                'displayRGCID', false, ...
                'pauseAfterEachRGCisRendered', ~true, ...
                'figureHandle', hFig, ...
                'videoOBJ', videoOBJ, ...
                'axesHandle', ax);
            pause
        end

    end % iCone

    if (visualizationParams.exportStagesAsPDF)
        stageName = 'Stage-1B';
        [hFig, ax] = setupFigure(stageName);
        [~, ~, XLims, YLims] = RGCRFconnector.visualizeWiringOfInputConeMosaicToRFcenters(...
            tmpRGCRFinputs, [], ...
            theInputConeMosaic, ...
            'titleString', stageName, ...
            'superimposeConeInputWiring', true, ...
            'displayRGCID', false, ...
            'pauseAfterEachRGCisRendered', ~true, ...
            'figureHandle', hFig, ...
            'videoOBJ', videoOBJ, ...
            'axesHandle', ax);

        NicePlot.exportFigToPDF(sprintf('%s-Stage1B.pdf', PDFfilename), hFig, 300);
    end

    
    % Find RGCs that end up having 0 cone inputs
    nonZeroInputRGCs = [];
    for iRGC = 1:numel(tmpRGCRFinputs)
        numberOfInputs = numel(tmpRGCRFinputs{iRGC});
        if (numberOfInputs > 0)
            nonZeroInputRGCs(numel(nonZeroInputRGCs)+1) = iRGC;
        end
    end
    % Keep a track of how many RGCs in the RGC lattice had 0 inputs.
    % These can be reused later
    availableZeroInputRGCsNum = numel(tmpRGCRFinputs)-numel(nonZeroInputRGCs);


    % For now just keep the non-zero cone input RGCs
    RGCRFinputsTmp = cell(1,numel(nonZeroInputRGCs));
    RGCRFweightsTmp = cell(1, numel(nonZeroInputRGCs));
    
    % Compute centroids and assigned weights to all cone inputs (all 1s).
    centroids = zeros(numel(nonZeroInputRGCs),2);
    parfor kkk = 1:numel(nonZeroInputRGCs)
        RGCRFinputsTmp{kkk} = tmpRGCRFinputs{nonZeroInputRGCs(kkk)};
        RGCRFweightsTmp{kkk} = RGCRFinputsTmp{kkk}*0+1;
        centroids(kkk,:) = RGCRFconnector.centroidsFromConeInputs({RGCRFinputsTmp{kkk}}, {RGCRFweightsTmp{kkk}}, allConePositions);
    end

    % Sort RFs according to distance from the center of the RGCmosaic
    mosaicCenter = mean(centroids,1);
    d = sum((bsxfun(@minus, centroids, mosaicCenter)).^2, 2);
    [~,idx] = sort(d);

    % Finalize the computed params
    RGCRFinputs = cell(1,numel(nonZeroInputRGCs));
    RGCRFweights = cell(1, numel(nonZeroInputRGCs));
    for kkk = 1:numel(idx)
        RGCRFinputs{kkk} = RGCRFinputsTmp{idx(kkk)};
        RGCRFweights{kkk} = RGCRFweightsTmp{idx(kkk)};
    end

    if ((visualizationParams.exportStagesAsPDF) || (visualizationParams.generateProgressionVideo))
        stageName = 'Stage-1';
        [hFig, ax] = setupFigure(stageName);
        [~, ~, XLims, YLims] = RGCRFconnector.visualizeWiringOfInputConeMosaicToRFcenters(...
            RGCRFinputs, RGCRFweights, ...
            theInputConeMosaic, ...
            'titleString', stageName, ...
            'superimposeConeInputWiring', true, ...
            'pauseAfterEachRGCisRendered', ~true, ...
            'figureHandle', hFig, ...
            'videoOBJ', videoOBJ, ...
            'axesHandle', ax);
         NicePlot.exportFigToPDF(sprintf('%s.pdf', PDFfilename), hFig, 300);

         if (visualizationParams.generateProgressionVideo)
            % Close video stream
            videoOBJ.close();
         end
    end

end


function [RGCRFpos, tmpRGCRFinputs, connectedConeIndices] = DoIt(...
        RGCRFpos, tmpRGCRFinputs, ...
        localConeToRGCDensityRatios, ...
        connectedConeIndices, ...
        allConePositions,  ...
        coneTypes, targetConeTypes)

    for iRGC = 1:size(RGCRFpos,1)

        % Find the localConeToRGCDensityRatios closest cones to each RGCRF
        [~,closestConeIndices] = pdist2(allConePositions, RGCRFpos(iRGC,:), ...
            '', 'smallest', floor(localConeToRGCDensityRatios(iRGC)));

        % Find which of these cones are of the accepted type and not
        % already connected
        idx = find(...
            (ismember(coneTypes(closestConeIndices), targetConeTypes)) & ...
            (~ismember(closestConeIndices, connectedConeIndices)));
     
        coneIndicesToConnectToThisRGC = closestConeIndices(idx);

        % Update the inputs and the position of this RGC
        tmpRGCRFinputs{iRGC} = coneIndicesToConnectToThisRGC;
        RGCRFpos(iRGC,:) = mean(allConePositions(coneIndicesToConnectToThisRGC,:),1);

        % Keep track of connected cone indices in this stage
        connectedConeIndices = cat(1, connectedConeIndices, coneIndicesToConnectToThisRGC);
    end

end


function[RGCRFpos, tmpRGCRFinputs, unconnectedConeIndices] = DoItOLD(...
        RGCRFpos, tmpRGCRFinputs, unconnectedConeIndices, ...
        closestMRGCIndices, allConePositions, sortedConeIndices, ...
        coneTypes, targetConeTypes)

    for iCone = 1:numel(sortedConeIndices)
        theConeIndex = sortedConeIndices(iCone);
        theConeRFPosition = allConePositions(theConeIndex,:);

        if (~ismember(coneTypes(theConeIndex), targetConeTypes)) 
            % Ignore cones other than the target cone type
            continue;
        end

        % Find the RGC that is closest to this cone
        theRGCindex = closestMRGCIndices(theConeIndex);

        % How may inputs does it have?
        inputConeIndices = tmpRGCRFinputs{theRGCindex};
        switch (numel(inputConeIndices))
            case 0
                % Bingo. Nearest RGC has 0 inputs. Connect to it.
                tmpRGCRFinputs{theRGCindex} = theConeIndex;
                % Update the position of this RGCRF
                RGCRFpos(theRGCindex,:) = theConeRFPosition;
                
            case 1
                % Nearest RGC has 1 input already. So the cone will remain
                % unconnected. Added to the list of unconnected cones
                unconnectedConeIndices(numel(unconnectedConeIndices)+1) = theConeIndex;
                continue;
            otherwise
                error('How can this RGC have %d input at this stage?', numel(inputConeIndices));
        end
    end
end





function [hFig, ax] = setupFigure(FigureTitle)
    hFig = figure(100); clf;
    set(hFig, 'Name', FigureTitle, 'Color', [1 1 1], 'Position', [10 10 800 800]);
    ax = subplot('Position', [0.1 0.1 0.8 0.8]);
end
