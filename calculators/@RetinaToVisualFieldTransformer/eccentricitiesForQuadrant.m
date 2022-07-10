% Retrieve valid eccentricities for a retinal quadrant
function [horizontalEccDegs, verticalEccDegs, eccDegsForPlotting] = ...
    eccentricitiesForQuadrant(retinalQuadrant, whichEye, maxEcc)

    assert(ismember(retinalQuadrant, RetinaToVisualFieldTransformer.validRetinalQuadrants), ...
        sprintf('Unknown retinal quadrant: ''%s''\n.', retinalQuadrant));

    switch (retinalQuadrant)

        case RetinaToVisualFieldTransformer.temporalRetinaQuadrant
            if (strcmp(whichEye, RetinaToVisualFieldTransformer.leftEye))
                horizontalEccDegs = 0:40;
            else
                horizontalEccDegs = -(0:40);
            end

            % Only analyze up to maxEcc
            idx = find(abs(horizontalEccDegs)<=maxEcc);
            horizontalEccDegs = horizontalEccDegs(idx);

            verticalEccDegs = zeros(1,numel(horizontalEccDegs));
            eccDegsForPlotting = abs(horizontalEccDegs);

        case RetinaToVisualFieldTransformer.nasalRetinaQuadrant
            if (strcmp(whichEye, RetinaToVisualFieldTransformer.leftEye))
                horizontalEccDegs = -(0:40);
            else
                horizontalEccDegs = 0:40;
            end

            % Only analyze up to maxEcc
            idx = find(abs(horizontalEccDegs)<=maxEcc);
            horizontalEccDegs = horizontalEccDegs(idx);

            verticalEccDegs = zeros(1,numel(horizontalEccDegs));
            eccDegsForPlotting = abs(horizontalEccDegs);

        case RetinaToVisualFieldTransformer.inferiorRetinaQuadrant
            verticalEccDegs = -(0:25);

            % Only analyze up to maxEcc
            idx = find(abs(verticalEccDegs)<=maxEcc);
            verticalEccDegs  = verticalEccDegs(idx);

            horizontalEccDegs = zeros(1,numel(verticalEccDegs));
            eccDegsForPlotting = abs(verticalEccDegs);

        case RetinaToVisualFieldTransformer.superiorRetinaQuadrant
            verticalEccDegs = (0:25);

            % Only analyze up to maxEcc
            idx = find(abs(verticalEccDegs)<=maxEcc);
            verticalEccDegs  = verticalEccDegs(idx);
            
            horizontalEccDegs = zeros(1,numel(verticalEccDegs));
            eccDegsForPlotting = abs(verticalEccDegs);
  
    end
end