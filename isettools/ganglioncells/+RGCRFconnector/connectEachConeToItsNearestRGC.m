function  [RGCRFinputs, RGCRFweights, availableZeroInputRGCsNum] = ...
    connectEachConeToItsNearestRGC(coneRFpos, coneRFspacings, coneTypes, ...
                                   RGCRFpos, chromaticSpatialVarianceTradeoff, ...
                                   sequentialConeTypeWiring, ...
                                   maxNearbyRGCsNum, visualizeProgression)
% Connect each of a set of input cones (L- or M-) to a set of RGC RF centers
%
% Syntax:
%   [RGCRFinputs, RGCRFweights, availableZeroInputRGCsNum] = ...
%        RGCRFconnector.connectEachConeToItsNearestRGC(...
%                    coneRFpos, coneRFspacings, coneTypes, ...
%                    RGCRFpos, chromaticSpatialVarianceTradeoff)
%
% Description:
%   Connect each of a set of input cones (L- or M-) to a set of RGC RF centers
%
% Inputs:
%    coneRFpos                  - [N x 2] matrix of (x,y) positions of N input cones
%    coneRFspacings             - [N x 1] vector of spacings of N input cones
%    coneTypes                  - [N x 1] vector of types of N input cones
%    RGCRFpos                   - [M x 2] matrix of (x,y) positions of M target RGC RF centers
%    chromaticSpatialVarianceTradeoff  - Chromatic-SpatialVariance tradefoff
%    sequentialConeTypeWiring   - Logical. Whether to wire cone types sequentially or not
%    maxNearbyRGCsNum           - Scalar. Max number of nearby RGCs to look for assignment, if
%                                 a cone cannot be assigned to its closest
%                                 RGC because it already has an input cone
%    visualizeProgression       - Boolean. Show the different steps
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

    % Sort cones according to distance from the center of the cone mosaic
    mosaicCenter = mean(coneRFpos,1);
    d = sum((bsxfun(@minus, coneRFpos, mosaicCenter)).^2, 2);
    [~,sortedConeIndices] = sort(d);

    % Find the closest RGC to each cone
    [~,closestMRGCIndices] = pdist2(RGCRFpos, coneRFpos, '', 'smallest', 1);

    % Store the cone inputs in a temporary array
    tmpRGCRFinputs = cell(1,size(RGCRFpos,1));
    
    % Pass 1. Only connect a cone to its nearest RGC if and only if that
    % RGC has 0 inputs. Also, keep a track of cone indices that were not 
    % connected because the nearest RGC was already connected to another cone
    unconnectedConeIndices = [];
    
    if (visualizeProgression)
        % Video setup
        if (sequentialConeTypeWiring)
            videoOBJ = VideoWriter(sprintf('WiringSteps_W%1.2f_SequentialConeTypeWiring', chromaticSpatialVarianceTradeoff), 'MPEG-4');
        else
            videoOBJ = VideoWriter(sprintf('WiringSteps_W%1.2f_NonSequentialConeTypeWiring', chromaticSpatialVarianceTradeoff), 'MPEG-4');
        end

        videoOBJ.FrameRate = 10;
        videoOBJ.Quality = 100;
        videoOBJ.open();

        hFig = figure(100); clf;
        set(hFig, 'Name', 'Stage 1', 'Color', [1 1 1]);
        ax = subplot('Position', [0.1 0.1 0.8 0.8]);
    end

    if (sequentialConeTypeWiring)
        % Connect the M-cones during the first pass
        targetConeTypes = [cMosaic.MCONE_ID];
        [RGCRFpos, tmpRGCRFinputs, unconnectedConeIndices] = DoIt(...
            RGCRFpos, tmpRGCRFinputs, unconnectedConeIndices, ...
            closestMRGCIndices, coneRFpos, sortedConeIndices, ...
            coneTypes, targetConeTypes);
    
        % Connect the L-cones during the second pass
        targetConeTypes = [cMosaic.LCONE_ID];
        [RGCRFpos, tmpRGCRFinputs, unconnectedConeIndices] = DoIt(...
            RGCRFpos, tmpRGCRFinputs, unconnectedConeIndices, ...
            closestMRGCIndices, coneRFpos, sortedConeIndices, ...
            coneTypes, targetConeTypes);
    else
        % Connect both L and M-cones at the same pass
        targetConeTypes = [cMosaic.LCONE_ID cMosaic.MCONE_ID];
        [RGCRFpos, tmpRGCRFinputs, unconnectedConeIndices] = DoIt(...
            RGCRFpos, tmpRGCRFinputs, unconnectedConeIndices, ...
            closestMRGCIndices, coneRFpos, sortedConeIndices, ...
            coneTypes, targetConeTypes);
    end


    if (visualizeProgression)
        stageName = 'Stage-1';
        [hFig, ax] = setupFigure(stageName);
        [~, ~, XLims, YLims] = RGCRFconnector.visualizeWiringOfInputConeMosaicToRFcenters(...
            tmpRGCRFinputs, [], ...
            coneRFpos, coneRFspacings, coneTypes, ...
            'titleString', stageName, ...
            'pauseAfterEachRGCisRendered', false, ...
            'superimposeConeInputWiring', true, ...
            'figureHandle', hFig, ...
            'videoOBJ', videoOBJ, ...
            'axesHandle', ax);
    end

    
    fprintf('%d out of %d L/M cones could not be connected to their closest RGC because the closest RGC already had one L/M cone attached\n', ...
        numel(unconnectedConeIndices), numel((coneTypes ~= cMosaic.SCONE_ID)));
    fprintf('Trying to connect them to nearby RGCs, minimizing the cost\n');
    


    % Pass 2. Deal with unconnected cone indices
    for iCone = 1:numel(unconnectedConeIndices)
        theConeIndex = unconnectedConeIndices(iCone);
        theConeRFPosition = coneRFpos(theConeIndex,:);

        % Find the maxNearbyRGCsNum neigboring RGCs 
        [~, nearestRGCindices] = pdist2(RGCRFpos, theConeRFPosition, '', 'smallest', maxNearbyRGCsNum);

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
    
            for iRGC = 1:numel(nearestRGCindices)
                % Retrieve this cell's cone inputs
                theRGCindex = nearestRGCindices(iRGC);
                inputConeIndices = tmpRGCRFinputs{theRGCindex};
    
                % To test the new cost, add this cone as an input
                inputConeIndices(numel(inputConeIndices)+1) = theConeIndex;
                inputConePositions = coneRFpos(inputConeIndices,:);
                inputConeSpacings = coneRFspacings(inputConeIndices);
                inputConeTypes = coneTypes(inputConeIndices);
                inputConeWeights = ones(1, numel(inputConeIndices));
    
                % Compute new cost, had this cone been assigned to this RGC
                projectedCosts(iRGC) = RGCRFconnector.costToMaintainInputCones(...
                    chromaticSpatialVarianceTradeoff, ...
                    inputConePositions, inputConeSpacings, ...
                    inputConeTypes, inputConeWeights);
            end % iRGC

            % Find the min projected cost
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
                            d = pdist2(coneRFpos(theInputConeIndicesForThisRGC,:), theConeRFPosition);
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
            RGCRFpos(targetRGCindex,:) = RGCRFconnector.centroidsFromConeInputs({tmpRGCRFinputs{targetRGCindex}}, [], coneRFpos);
        end

        if (visualizeProgression)
            stageName = 'Stage-2';
            [hFig, ax] = setupFigure(stageName);
            [~, ~, XLims, YLims] = RGCRFconnector.visualizeWiringOfInputConeMosaicToRFcenters(...
                tmpRGCRFinputs, [], ...
                coneRFpos, coneRFspacings, coneTypes, ...
                'titleString', stageName, ...
                'superimposeConeInputWiring', true, ...
                'displayRGCID', false, ...
                'pauseAfterEachRGCisRendered', ~true, ...
                'figureHandle', hFig, ...
                'videoOBJ', videoOBJ, ...
                'axesHandle', ax);
        end

    end % iCone

    
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
        centroids(kkk,:) = RGCRFconnector.centroidsFromConeInputs({RGCRFinputsTmp{kkk}}, {RGCRFweightsTmp{kkk}}, coneRFpos);
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

    if (visualizeProgression)
        stageName = 'Stage-3';
        [hFig, ax] = setupFigure(stageName);
        [~, ~, XLims, YLims] = RGCRFconnector.visualizeWiringOfInputConeMosaicToRFcenters(...
            RGCRFinputs, RGCRFweights, ...
            coneRFpos, coneRFspacings, coneTypes, ...
            'titleString', stageName, ...
            'superimposeConeInputWiring', ~true, ...
            'pauseAfterEachRGCisRendered', ~true, ...
            'figureHandle', hFig, ...
            'videoOBJ', videoOBJ, ...
            'axesHandle', ax);
    end
   
    videoOBJ.close();
end



function[RGCRFpos, tmpRGCRFinputs, unconnectedConeIndices] = DoIt(...
        RGCRFpos, tmpRGCRFinputs, unconnectedConeIndices, ...
        closestMRGCIndices, coneRFpos, sortedConeIndices, ...
        coneTypes, targetConeTypes)

    for iCone = 1:numel(sortedConeIndices)
        theConeIndex = sortedConeIndices(iCone);
        theConeRFPosition = coneRFpos(theConeIndex,:);

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
