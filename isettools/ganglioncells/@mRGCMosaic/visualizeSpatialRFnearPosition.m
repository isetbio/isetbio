function visualizeSpatialRFnearPosition(obj, ...
    targetRGCposition, targetCenterConesNum, ...
    targetCenterConeMajorityType, varargin)

    % Parse optional input
    p = inputParser;
    p.addParameter('pdfFileName', '', @ischar);
    p.parse(varargin{:});
    pdfFileName = p.Results.pdfFileName;

    % Customize PDF filename
    if (~isempty(pdfFileName))
        if (isnan(targetCenterConeMajorityType))
            pdfPostFix = sprintf('AtPosition_%2.2f_%2.2f_CenterConesNum_%d_mixedLM', ...
                    targetRGCposition(1), targetRGCposition(2), targetCenterConesNum);
        else
            switch (targetCenterConeMajorityType)
                case cMosaic.LCONE_ID
                    pdfPostFix = sprintf('AtPosition_%2.2f_%2.2f_CenterConesNum_%d_LconeDominated.pdf', ...
                        targetRGCposition(1), targetRGCposition(2), targetCenterConesNum);
                case cMosaic.MCONE_ID
                    pdfPostFix = sprintf('AtPosition_%2.2f_%2.2f_CenterConesNum_%d_MconeDominated.pdf', ...
                        targetRGCposition(1), targetRGCposition(2), targetCenterConesNum);
                    
            end
        end
    
        pdfFileName = strrep(pdfFileName, '.pdf', pdfPostFix);
    end

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
        fprintf('Found mosaic with %d center cones near the target position (%2.2f,%2.2f) at %2.2f,%2.2f\n', ...
            targetCenterConesNum, ...
            targetRGCposition(1), targetRGCposition(2), ...
            obj.rgcRFpositionsDegs(theCurrentRGCindex,1), obj.rgcRFpositionsDegs(theCurrentRGCindex,2));
    end

    obj.visualizeSpatialRFofRGCwithIndex(theCurrentRGCindex, pdfFileName);
end