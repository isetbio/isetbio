function theVisualizedRGCindex = visualizeRetinalConePoolingRFmapNearPosition(obj, ...
    targetRGCposition, targetCenterConesNum, ...
    targetCenterConeMajorityType,  ...
    varargin)

    % Parse optional input
    p = inputParser;
    p.addParameter('tickSeparationArcMin', 6, @isscalar);
    p.addParameter('normalizedPeakSurroundSensitivity', 0.4, @isscalar);
    p.addParameter('withFigureFormat', [], @(x)(isempty(x)||(isstruct(x))));
    p.addParameter('theAxes', [], @(x)(isempty(x)||(iscell(x))));
    p.addParameter('reverseXDir', false, @islogical);
    p.addParameter('gridlessLineWeightingFunctions', false, @islogical);

    p.parse(varargin{:});
    
    ff = p.Results.withFigureFormat;
    tickSeparationArcMin = p.Results.tickSeparationArcMin;
    normalizedPeakSurroundSensitivity = p.Results.normalizedPeakSurroundSensitivity;
    theAxes = p.Results.theAxes;
    reverseXDir = p.Results.reverseXDir;
    gridlessLineWeightingFunctions = p.Results.gridlessLineWeightingFunctions;

    if (isempty(targetRGCposition))
        theVisualizedRGCindex = input('Enter index of target RGC : ');
    else
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
        end
    end

    obj.visualizeRetinalConePoolingRFmapOfRGCwithIndex(theVisualizedRGCindex, ...
            'theAxes', theAxes, ...
            'tickSeparationArcMin', tickSeparationArcMin, ...
            'normalizedPeakSurroundSensitivity', normalizedPeakSurroundSensitivity, ...
            'withFigureFormat', ff, ...
            'reverseXDir', reverseXDir, ...
            'gridlessLineWeightingFunctions', gridlessLineWeightingFunctions);
           
end