function divergeConeOutputsToMultipleNearbyRGCs(obj, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('RcToRGCseparationRatio', [], @(x)((isempty(x))||(isscalar(x)&&(x>=1))));
    p.parse(varargin{:});
    
    % Remove any negative weights (which indicate overlapping cone weights)
    obj.coneConnectivityMatrix(obj.coneConnectivityMatrix<0) = 0;

    if (~isempty(p.Results.RcToRGCseparationRatio))
        obj.wiringParams.RcToRGCseparationRatio = p.Results.RcToRGCseparationRatio;
    end

    if (obj.wiringParams.RcToRGCseparationRatio == 1.0)
        return;
    end

    rgcsNum = size(obj.coneConnectivityMatrix,2);
    conesNum = size(obj.coneConnectivityMatrix,1);
    
    allConePositions = obj.inputConeMosaic.coneRFpositionsMicrons;
    
    activatedConeIndices = [obj.inputConeMosaic.mConeIndices(:); obj.inputConeMosaic.lConeIndices(:)];
    theFullConnectivityMatrix = cell(rgcsNum,2);
 
    for iRGC = 1:rgcsNum
        % The centroid of the RGC
        theRGCCentroid = obj.RGCRFcentroidsFromInputs(iRGC,:);
        
        % Indices of non-overlapping input cones
        inputConeIndices = find(squeeze(obj.coneConnectivityMatrix(:, iRGC))>0);
        % Weights of non-overlapping input cones
        inputConeWeights = full(obj.coneConnectivityMatrix(inputConeIndices, iRGC));
     
        % The pooling radius of the RGC
        % meanInputConeSpacing = mean(obj.inputConeMosaic.coneRFspacingsMicrons(inputConeIndices));
        
        overlapCharacteristicRadius = 0.204*sqrt(2.0)*obj.localRGCRFspacingsMicrons(iRGC) * (obj.wiringParams.RcToRGCseparationRatio);
        overlapRadius = 0.7*obj.localRGCRFspacingsMicrons(iRGC) * obj.wiringParams.RcToRGCseparationRatio;
        
        % Find cone indices within pooling radius distance from centroid
        distances = RGCconnector.pdist2(allConePositions(activatedConeIndices,:), theRGCCentroid);
        
        % Find the indices of the overlapping cones
        idx = find(distances <= overlapRadius);
        coneIndicesWithinOverlapRadius = activatedConeIndices(idx);
        coneDistancesWithinOverlapRadius = distances(idx);
        coneWeightsWithinOverlapRadius = exp(-(coneDistancesWithinOverlapRadius/overlapCharacteristicRadius).^2);

        % Find which of the coneIndicesWithinPoolingRadius are the main (non-overlapping cones)
        [isMainCone, ia] = ismember(coneIndicesWithinOverlapRadius, inputConeIndices);
        
        % Find the max Gaussian weight that would have been assigned to the non-overlapping  cones
        [overlappingConesFound,ib] = ismember(inputConeIndices, coneIndicesWithinOverlapRadius);

         
        if (~any(overlappingConesFound))
            % No cones exist that are overlapping and non main ones
            theFullConnectivityMatrix{iRGC}{1} = inputConeIndices;
            theFullConnectivityMatrix{iRGC}{2} = inputConeWeights ;
            continue;
        end

        scalingFactor = max(coneWeightsWithinOverlapRadius(ib(ib>0)));
        
        for iCone = 1:numel(coneIndicesWithinOverlapRadius)
            if (isMainCone(iCone))
                % Non-overlapping cone, so keep original weight
                coneWeightsWithinOverlapRadius(iCone) = inputConeWeights(ia(iCone));
            else
                % Overlapping cone. Set the weights to negative polarity so as to
                % discriminate between main and overlapping cones
                w = min([1 coneWeightsWithinOverlapRadius(iCone)/scalingFactor]);
                if (w > 1) || ( w < 0)
                    error('How can this be?')
                end
                coneWeightsWithinOverlapRadius(iCone) = -w;
            end
        end

        % Modify connectivity matrix to encode the weights of the overlaping cone inputs. 
        theFullConnectivityMatrix{iRGC}{1} = coneIndicesWithinOverlapRadius;
        theFullConnectivityMatrix{iRGC}{2} = coneWeightsWithinOverlapRadius;
    end
    
    % Finalize connectivity matrix
    for iRGC = 1:rgcsNum
        inputConeIndices = theFullConnectivityMatrix{iRGC}{1};
        inputConeWeights = theFullConnectivityMatrix{iRGC}{2};
        if (sum(inputConeWeights(:)) == 0)
            error('How can this be?')
        end
        obj.coneConnectivityMatrix(inputConeIndices, iRGC) = inputConeWeights;
    end

    
    % Find indices of cones that are connected to RGCs
    connectedConeIndices = find(sum(abs(obj.coneConnectivityMatrix),2)>0);
    
    % Sum of the output of each cone to all RGCs must equal 1.0;
    totalOutputForEachConnectedCone = sum(abs(obj.coneConnectivityMatrix(connectedConeIndices,:)),2);
    obj.coneConnectivityMatrix(connectedConeIndices,:) = ...
        bsxfun(@times, ...
        obj.coneConnectivityMatrix(connectedConeIndices,:), ...
        1./totalOutputForEachConnectedCone);
    
    fprintf('There are %d cones (%d of which are connected to %d RGCs\n', conesNum, numel(connectedConeIndices), rgcsNum);
    totalInputForEachRGC = sum(abs(obj.coneConnectivityMatrix(connectedConeIndices,:)),1);
    totalOutputForEachConnectedCone = sum(abs(obj.coneConnectivityMatrix(connectedConeIndices,:)),2);
    fprintf('max input across all %d RGCs: %f\n', numel(totalInputForEachRGC), max(full(totalInputForEachRGC)));
    fprintf('min input across all %d RGCs: %f\n', numel(totalInputForEachRGC), min(full(totalInputForEachRGC)));
    fprintf('max output across all %d connected cones: %f\n', numel(totalOutputForEachConnectedCone), max(full(totalOutputForEachConnectedCone)));
    fprintf('min output across all %d connected cones: %f\n', numel(totalOutputForEachConnectedCone), min(full(totalOutputForEachConnectedCone)));
end

