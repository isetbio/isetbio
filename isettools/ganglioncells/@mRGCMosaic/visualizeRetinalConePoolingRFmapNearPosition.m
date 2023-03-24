function theVisualizedRGCindex = visualizeRetinalConePoolingRFmapNearPosition(obj, ...
    targetRGCposition, targetCenterConesNum, ...
    targetCenterConeMajorityType, theAxes, ...
    varargin)

     % Parse optional input
     p = inputParser;
     p.addParameter('withFigureFormat', [], @(x)(isempty(x)||(isstruct(x))));
     p.parse(varargin{:});
     ff = p.Results.withFigureFormat;


    % Find the RGC index that best matches the target criteria
    [targetCenterConesNumNotMatched, theVisualizedRGCindex] = obj.indexOfRGCNearPosition(...
        targetRGCposition, targetCenterConesNum, targetCenterConeMajorityType);

    if (targetCenterConesNumNotMatched)
        fprintf(2, 'Could not find an RGC with %d center cones\n', targetCenterConesNum);
        return;
    else
        fprintf('Found mosaic with %d center cones near the target position (%2.2f,%2.2f) at %2.2f,%2.2f\n', ...
            targetCenterConesNum, ...
            targetRGCposition(1), targetRGCposition(2), ...
            obj.rgcRFpositionsDegs(theVisualizedRGCindex,1), obj.rgcRFpositionsDegs(theVisualizedRGCindex,2));
    
        obj.visualizeRetinalConePoolingRFmapOfRGCwithIndex(theVisualizedRGCindex, theAxes, ...
            'withFigureFormat', ff);
    end
end