function [sortedIndices, sortedEccs] = sortDestinationRFsBasedOnOptimizationCenter(obj,unsortedIndices, varargin)

    p = inputParser;
    p.addParameter('ignoreInfCentroids', false, @islogical);
    p.addParameter('averageEcc', false, @islogical);
    p.parse(varargin{:});
    ignoreInfCentroids = p.Results.ignoreInfCentroids;
    averageEcc = p.Results.averageEcc;

    switch (obj.wiringParams.optimizationCenter)
        case 'origin'
            ecc = sqrt(sum(obj.destinationRFcentroidsFromInputs.^2,2));

        case 'latticeCenter'
            if (isempty(obj.sourceLatticeCenter))
                obj.sourceLatticeCenter = mean(obj.sourceLattice.RFpositionsMicrons,1);
            end
            diff = bsxfun(@minus, obj.destinationRFcentroidsFromInputs, obj.sourceLatticeCenter);
            ecc = sqrt(sum(diff.^2,2));

        case 'localSpacing'
            ecc = obj.destinationRFspacingsFromCentroids; 

        case 'localConeToRGCdensityRatio'
            % Compute the source-to-destination density ratio map at the current
            % RFpos of the destination lattice
            ecc = obj.sourceToDestinationDensityRatioMap();

        otherwise
            error('Unknown wiringParams.optimizationCenter: ''%s''.', obj.wiringParams.optimizationCenter);
    end % switch

    
    if (averageEcc)
        eccAllDestinationRFs = ecc;
        ecc = zeros(1,numel(unsortedIndices));
        parfor iDestinationRF = 1:numel(unsortedIndices)
            destinationRF = unsortedIndices(iDestinationRF);
            nearbyDestinationRFindicesForDestinationRF = obj.indicesOfNeighboringDestinationRFs(...
                destinationRF, ...
                'ignoreInfCentroids', ignoreInfCentroids);
            nearbyDestinationRFindicesForDestinationRF(numel(nearbyDestinationRFindicesForDestinationRF)+1) = destinationRF;
            nearbyEccs = eccAllDestinationRFs(nearbyDestinationRFindicesForDestinationRF);
            nearbyEccs = nearbyEccs(:);
            idx = find(~isinf(nearbyEccs));
            ecc(iDestinationRF) = median(nearbyEccs(idx));
        end
    else
        ecc = ecc(unsortedIndices);
    end

   % Compute sorted indices of destination RFs in increasing eccentricity
   % or local spacing
   [sortedEccs, sortedIndices] = sort(ecc, 'ascend');
end
