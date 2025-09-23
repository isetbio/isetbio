function s = singleCellConnectivityStats(obj, theRGCindex, theSubregion, varargin)
    
    % Parse input
    p = inputParser;
    p.addParameter('inputConeIndicesOnly', false, @islogical);
    p.addParameter('minConeWeightIncluded', mRGCMosaic.sensitivityAtPointOfOverlap, @isscalar);
    p.addParameter('coneTypesConsidered', [cMosaic.LCONE_ID cMosaic.MCONE_ID]);
    p.addParameter('warnIfCenterConeInputNumerosityDiffersFromExclusiveOne', true, @islogical);

    p.parse(varargin{:});
    inputConeIndicesOnly = p.Results.inputConeIndicesOnly;
    minConeWeightIncluded = p.Results.minConeWeightIncluded;
    coneTypesConsidered = p.Results.coneTypesConsidered;
    warnIfCenterConeInputNumerosityDiffersFromExclusiveOne = p.Results.warnIfCenterConeInputNumerosityDiffersFromExclusiveOne;

    assert(ismember(theSubregion, {'center', 'surround', 'surround-center'}), 'theSubregion must be either ''center'' or ''surround'' or ''surround-center''.');

    % Retrieve this cell's # of center cone indices
    if (strcmp(theSubregion, 'center'))
        connectivityVector = full(squeeze(obj.rgcRFcenterConeConnectivityMatrix(:, theRGCindex)));
        theInputConeIndices = find(...
            (connectivityVector >=  minConeWeightIncluded) & ...
            (ismember(obj.inputConeMosaic.coneTypes, coneTypesConsidered)));


        if (~isempty(obj.exclusivelyConnectedInputConeIndicesNum))&&(numel(theInputConeIndices) ~= obj.exclusivelyConnectedInputConeIndicesNum(theRGCindex))
            if (warnIfCenterConeInputNumerosityDiffersFromExclusiveOne)
                fprintf('**Different number of input cones at minConeWeightIncluded = %f ((overlapping:%d vs exclusive:%d)\n', minConeWeightIncluded, numel(theInputConeIndices), obj.exclusivelyConnectedInputConeIndicesNum(theRGCindex));
            end
        end

        if (numel(theInputConeIndices) == 0)
            minConeWeightIncluded
            [numel(theInputConeIndices) obj.exclusivelyConnectedInputConeIndicesNum(theRGCindex)]
                theRGCindex
                pause
                obj.rgcRFpositionsDegs(theRGCindex,:)
                numel(find(connectivityVector>0))
                [min(connectivityVector) max(connectivityVector)]
                numel(theInputConeIndices)
                theInputConeIndices
                full(connectivityVector(theInputConeIndices, 1))
                 max(connectivityVector) 
        end

    else
        if (strcmp(theSubregion, 'surround'))
            % For the surround, we include weights down to the
            % minConeWeightIncluded scaled by the surround peak amplitude
            connectivityVector = full(squeeze(obj.rgcRFsurroundConeConnectivityMatrix(:, theRGCindex)));
            theInputConeIndices = find(connectivityVector > max(connectivityVector)*minConeWeightIncluded);

        elseif (strcmp(theSubregion, 'surround-center'))
            % Here we return the surround cone indices for which the composite surround-center weights > 0
            connectivityVector = full(squeeze(obj.rgcRFcenterConeConnectivityMatrix(:, theRGCindex)));
            allCenterConnectedConeIndices = find(connectivityVector>0.0001);
            connectivityVector = full(squeeze(obj.rgcRFsurroundConeConnectivityMatrix(:, theRGCindex)));
            allSurroundConnectedConeIndices = find(connectivityVector>0.0001);

            % cones connected to both center and surround with composite surround-center weight > 0
            theSurroundCenterCommonConeIndices = intersect(allCenterConnectedConeIndices, allSurroundConnectedConeIndices);
            theSurroundCenterCompositeConeWeights = ...
                obj.rgcRFsurroundConeConnectivityMatrix(theSurroundCenterCommonConeIndices, theRGCindex) - ...
                obj.rgcRFcenterConeConnectivityMatrix(theSurroundCenterCommonConeIndices, theRGCindex);

            % cones whose surround cone pooling weight exceeds that of the center
            theSurroundConeIndicesWithWeightsStrongerThanTheirCenterWeights = theSurroundCenterCommonConeIndices(theSurroundCenterCompositeConeWeights>0);

            % Cones exclusively connected to the surround
            theSurroundExclusiveConeIndices = allSurroundConnectedConeIndices(~ismember(allSurroundConnectedConeIndices, allCenterConnectedConeIndices));

            % Return the indices for which the surround-center composite weights are > 0
            theInputConeIndices = union(theSurroundConeIndicesWithWeightsStrongerThanTheirCenterWeights, theSurroundExclusiveConeIndices);
        end
    end
    
    if (inputConeIndicesOnly)
        s = theInputConeIndices;
        return;
    end

    if (isempty(theInputConeIndices))
        fprintf(2, 'cell #%d has zero %s inputs with weight greater than the specified min cone weight: %f, Returning [].\n', theRGCindex, theSubregion, minConeWeightIncluded)
        s = [];
        return;
    end

    theInputConeWeights = connectivityVector(theInputConeIndices);
    theInputConeTypes = obj.inputConeMosaic.coneTypes(theInputConeIndices);

    s.inputConesNum = numel(theInputConeIndices);
    s.netWeights(cMosaic.LCONE_ID) = sum(theInputConeWeights((theInputConeTypes == cMosaic.LCONE_ID)));
    s.netWeights(cMosaic.MCONE_ID) = sum(theInputConeWeights((theInputConeTypes == cMosaic.MCONE_ID)));
    s.netWeights(cMosaic.SCONE_ID) = sum(theInputConeWeights((theInputConeTypes == cMosaic.SCONE_ID)));
    s.netWeights = s.netWeights / sum(s.netWeights(:));

    [~, dominantConeType] = max(s.netWeights);
    switch (dominantConeType)
        case 1
            s.dominantConeType = cMosaic.LCONE_ID;
        case 2
            s.dominantConeType = cMosaic.MCONE_ID;
        case 3
            s.dominantConeType = cMosaic.SCONE_ID;
    end
    
end