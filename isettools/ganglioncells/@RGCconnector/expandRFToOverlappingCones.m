function expandRFToOverlappingCones(obj, varargin)

    % Parse input
    %p = inputParser;
    %p.addParameter('generateProgressVideo', false, @islogical);
    %p.addParameter('optimizationCenter', 'patchCenter', @(x)(ismember(x, {'patchCenter', 'visualFieldCenter'})));

    %p.parse(varargin{:});
    %generateProgressVideo = p.Results.generateProgressVideo;
    %optimizationCenter = p.Results.optimizationCenter;
    
    if (obj.wiringParams.rfOverlapFactor <= 0)
        return;
    end
    
    rgcsNum = size(obj.RGCRFcentroidsFromInputs,1);
    allConePositions =  obj.inputConeMosaic.coneRFpositionsMicrons;
    
    activatedConeIndices = [obj.inputConeMosaic.mConeIndices(:); obj.inputConeMosaic.lConeIndices(:)];
    theFullConnectivityMatrix = cell(rgcsNum,2);
    
   
    for iRGC = 1:rgcsNum
        % The centroid of the RGC
        theRGCCentroid = obj.RGCRFcentroidsFromInputs(iRGC,:);
        
        % The pooling radius of the RGC
        poolingSigma = obj.localRGCRFspacingsMicrons(iRGC)*obj.wiringParams.rfOverlapFactor;
        poolingRadius = poolingSigma * 3.0;
        
        % Find cone indices within pooling radius distance from centroid
        distances = RGCconnector.pdist2(allConePositions(activatedConeIndices,:), theRGCCentroid);
        
        % Find the indices of the overlapping cones
        idx = find(distances < poolingRadius);
        coneIndicesWithinPoolingRadius = activatedConeIndices(idx);
        
        coneDistancesWithinPoolingRadius = distances(idx);
        coneWeightsWithinPoolingRadius = exp(-0.5*(coneDistancesWithinPoolingRadius/poolingSigma).^2);
        
        % Indices of non-overlapping input cones
        connectedConeIndices = find(squeeze(obj.coneConnectivityMatrix(:, iRGC))>0);
        % Weights of non-overlapping input cones
        inputConeWeights = full(obj.coneConnectivityMatrix(connectedConeIndices, iRGC));
     
        % Find which of the coneIndicesWithinPoolingRadius are also the
        % main (non-overlapping cones)
        [isMainCone, ia] = ismember(coneIndicesWithinPoolingRadius, connectedConeIndices);
        
        % Find the max Gaussian weight that would have been assigned to the non-overlapping  cones
        [~,ib] = ismember(connectedConeIndices, coneIndicesWithinPoolingRadius);
        scalingFactor = max(coneWeightsWithinPoolingRadius(ib));
        
        for iCone = 1:numel(coneIndicesWithinPoolingRadius)
            if (isMainCone(iCone))
                % Non-overlapping cone, so keep original weight
                coneWeightsWithinPoolingRadius(iCone) = inputConeWeights(ia(iCone));
            else
                % Overlapping cone. Set the weights to negative polarity so as to
                % discriminate between main and overlapping cones
                w = coneWeightsWithinPoolingRadius(iCone)/scalingFactor;
                if (w > 1) || ( w < 0)
                    error('How can this be?')
                end
                coneWeightsWithinPoolingRadius(iCone) = -w;
            end
        end
        
        % Modify connectivity matrix to encode the weights of the overlaping cone inputs. 
        theFullConnectivityMatrix{iRGC}{1} = coneIndicesWithinPoolingRadius;
        theFullConnectivityMatrix{iRGC}{2} = coneWeightsWithinPoolingRadius;
    end
    
    for iRGC = 1:rgcsNum
        connectedConeIndices = theFullConnectivityMatrix{iRGC}{1};
        connectedConeWeights = theFullConnectivityMatrix{iRGC}{2};
        obj.coneConnectivityMatrix(connectedConeIndices, iRGC) = connectedConeWeights;
    end
    
    
    
end

