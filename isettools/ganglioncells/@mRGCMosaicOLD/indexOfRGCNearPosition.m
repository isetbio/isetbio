function [targetCenterConesNumNotMatched, theCurrentRGCindex] = indexOfRGCNearPosition(obj, ...
    targetRGCposition, targetCenterConesNum, targetCenterConeMajorityType)

    targetCenterConesNumNotMatched = true;
    theCurrentRGCindex = [];
    iRGC = 0;

    beVerbose = false;

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
                if (beVerbose)
                fprintf('Nope: cell at distance %2.3f has %d cones in the center with majority type %d\n', ...
                    sortedDistances(iRGC), sum(theCenterConeTypeNum), theMajorityConeType);
                end
            else
                targetCenterConesNumNotMatched = false;
            end
        else
            if (sum(theCenterConeTypeNum) ~= targetCenterConesNum) || (theMajorityConeType ~= targetCenterConeMajorityType)
                if (beVerbose)
                fprintf('Nope: cell at distance %2.3f has %d cones in the center with majority type %d\n', ...
                    sortedDistances(iRGC), sum(theCenterConeTypeNum), theMajorityConeType);
                end
            else
                targetCenterConesNumNotMatched = false;
            end
        end
    end

    if (targetCenterConesNumNotMatched)
        fprintf(2, 'Did not find an RGC near (%2.1f, %2.1f) with %d center cones of type %d\n', ...
            targetRGCposition(1), targetRGCposition(2), targetCenterConesNum, targetCenterConeMajorityType)
    end

end
