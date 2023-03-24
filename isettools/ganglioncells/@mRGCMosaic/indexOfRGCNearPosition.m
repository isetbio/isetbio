function [targetCenterConesNumNotMatched, theCurrentRGCindex] = indexOfRGCNearPosition(obj, ...
    targetRGCposition, targetCenterConesNum, targetCenterConeMajorityType)

    targetCenterConesNumNotMatched = true;
    theCurrentRGCindex = [];
    iRGC = 0;

    d = (bsxfun(@minus, obj.rgcRFpositionsDegs, targetRGCposition)).^2;
    d = sqrt(sum(d,2));
    [sortedDistances,sortedRGCindices] = sort(d,'ascend');

    while ((targetCenterConesNumNotMatched) && (iRGC < numel(sortedRGCindices)))
        iRGC = iRGC+1;
        theCurrentRGCindex = sortedRGCindices(iRGC);

        [theCenterConeTypeWeights, theCenterConeTypeNum, theMajorityConeType, theCenterConeTypes] = ...
            obj.centerConeTypeWeights(theCurrentRGCindex);

        if (isnan(targetCenterConeMajorityType))
            if (sum(theCenterConeTypeNum) ~= targetCenterConesNum) || (~isnan(theMajorityConeType))
                fprintf('Nope: cell at distance %2.3f has %d cones in the center with majority type %d\n', ...
                    sortedDistances(iRGC), sum(theCenterConeTypeNum), theMajorityConeType);
            else
                targetCenterConesNumNotMatched = false;
            end
        else
            if (sum(theCenterConeTypeNum) ~= targetCenterConesNum) || (theMajorityConeType ~= targetCenterConeMajorityType)
                fprintf('Nope: cell at distance %2.3f has %d cones in the center with majority type %d\n', ...
                    sortedDistances(iRGC), sum(theCenterConeTypeNum), theMajorityConeType);
            else
                targetCenterConesNumNotMatched = false;
            end
        end
    end
end
