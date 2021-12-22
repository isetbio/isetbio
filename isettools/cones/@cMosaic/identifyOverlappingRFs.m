function [rfsToKeep, rfsToBeEliminated, overlapingRFindex] = identifyOverlappingRFs(xPos, yPos, ...
             RFpositionsMicrons, RFspacingsMicrons, maxSeparationForDeclaringOverlap)
    
    rfsNum  = size(RFpositionsMicrons,1);
    overlapingRFindex = zeros(1,rfsNum);

    parfor iRF = 1:rfsNum-1
        otherIndices = iRF+1:rfsNum;
        dd = sqrt(sum((RFpositionsMicrons(iRF,:)-RFpositionsMicrons(otherIndices,:)).^2,2));
        [minDistance, idx] = min(dd);
        otherRF = otherIndices(idx);
        meanSpacing = max([RFspacingsMicrons(iRF) RFspacingsMicrons(otherRF)]);
        if (minDistance < meanSpacing*maxSeparationForDeclaringOverlap)
            overlapingRFindex(iRF) = otherRF;
        end
    end % iRF

    rfsToBeEliminated = find(overlapingRFindex > 0);
    problematicNodesNum = numel(rfsToBeEliminated);
    rfsToKeep = setdiff(1:rfsNum, rfsToBeEliminated);
end
