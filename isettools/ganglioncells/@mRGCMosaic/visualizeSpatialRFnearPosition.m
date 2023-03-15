function theVisualizedRGCindex = visualizeSpatialRFnearPosition(obj, ...
    targetRGCposition, targetCenterConesNum, ...
    targetCenterConeMajorityType, theAxes, varargin)

     % Parse optional input
     p = inputParser;
     p.addParameter('withFigureFormat', [], @(x)(isempty(x)||(isstruct(x))));
     p.parse(varargin{:});
     ff = p.Results.withFigureFormat;


    d = (bsxfun(@minus, obj.rgcRFpositionsDegs, targetRGCposition)).^2;
    d = sqrt(sum(d,2));
    [sortedDistances,sortedRGCindices] = sort(d,'ascend');

    targetCenterConesNumNotMatched = true;
    iRGC = 0;
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

    if (targetCenterConesNumNotMatched)
        fprintf(2, 'Could not find an RGC with %d center cones\n', targetCenterConesNum)
        return;
    else
        theVisualizedRGCindex = theCurrentRGCindex;
        fprintf('Found mosaic with %d center cones near the target position (%2.2f,%2.2f) at %2.2f,%2.2f\n', ...
            targetCenterConesNum, ...
            targetRGCposition(1), targetRGCposition(2), ...
            obj.rgcRFpositionsDegs(theVisualizedRGCindex,1), obj.rgcRFpositionsDegs(theVisualizedRGCindex,2));
    end

    obj.visualizeSpatialRFofRGCwithIndex(theVisualizedRGCindex, theAxes, ...
        'withFigureFormat', ff);
end