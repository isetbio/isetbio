function expandRFsToOverlappingCones(obj, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('rfOverlapFactor', [], @(x)((isempty(x))||(isscalar(x)&&(x>=0)&&(x<=1))));
    p.parse(varargin{:});
    
    % Remove any negative weights (which indicate overlapping cone weights)
    obj.coneConnectivityMatrix(obj.coneConnectivityMatrix<0) = 0;

    if (~isempty(p.Results.rfOverlapFactor))
        obj.wiringParams.rfOverlapFactor = p.Results.rfOverlapFactor;
    end

    if (obj.wiringParams.rfOverlapFactor <= 0)
        return;
    end

    rgcsNum = size(obj.RGCRFcentroidsFromInputs,1);
    allConePositions = obj.inputConeMosaic.coneRFpositionsMicrons;
    
    activatedConeIndices = [obj.inputConeMosaic.mConeIndices(:); obj.inputConeMosaic.lConeIndices(:)];
    theFullConnectivityMatrix = cell(rgcsNum,2);
 
    for iRGC = 1:rgcsNum
        % The centroid of the RGC
        theRGCCentroid = obj.RGCRFcentroidsFromInputs(iRGC,:);
        
        % Indices of non-overlapping input cones
        connectedConeIndices = find(squeeze(obj.coneConnectivityMatrix(:, iRGC))>0);
        % Weights of non-overlapping input cones
        inputConeWeights = full(obj.coneConnectivityMatrix(connectedConeIndices, iRGC));
     
        % The pooling radius of the RGC
        meanInputConeSpacing = mean(obj.inputConeMosaic.coneRFspacingsMicrons(connectedConeIndices));
        overlapSigma = meanInputConeSpacing * obj.wiringParams.rfOverlapFactor;
        overlapRadius = overlapSigma * 3.0;
        
        % Find cone indices within pooling radius distance from centroid
        distances = RGCconnector.pdist2(allConePositions(activatedConeIndices,:), theRGCCentroid);
        
        % Find the indices of the overlapping cones
        idx = find(distances < overlapRadius);
        coneIndicesWithinOverlapRadius = activatedConeIndices(idx);
        coneDistancesWithinOverlapRadius = distances(idx);
        coneWeightsWithinOverlapRadius = exp(-0.5*(coneDistancesWithinOverlapRadius/overlapSigma).^2);

        % Find which of the coneIndicesWithinPoolingRadius are the main (non-overlapping cones)
        [isMainCone, ia] = ismember(coneIndicesWithinOverlapRadius, connectedConeIndices);
        
        % Find the max Gaussian weight that would have been assigned to the non-overlapping  cones
        [overlappingConesFound,ib] = ismember(connectedConeIndices, coneIndicesWithinOverlapRadius);

         
        if (~any(overlappingConesFound))
            % No cones exist that are overlapping and non main ones
            theFullConnectivityMatrix{iRGC}{1} = connectedConeIndices;
            theFullConnectivityMatrix{iRGC}{2} = inputConeWeights ;
            continue;
        end

        scalingFactor = max(coneWeightsWithinOverlapRadius(ib));
        
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
        connectedConeIndices = theFullConnectivityMatrix{iRGC}{1};
        connectedConeWeights = theFullConnectivityMatrix{iRGC}{2};
        obj.coneConnectivityMatrix(connectedConeIndices, iRGC) = connectedConeWeights;
    end


end

