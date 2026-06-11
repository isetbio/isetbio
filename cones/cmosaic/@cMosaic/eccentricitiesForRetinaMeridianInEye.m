% Static method to return signed horizontal and vertical
% eccentricities corresponding to radial eccentricities specified
% on one of the 4 principal retinal meridians and eye
function [horizontalEcc, verticalEcc] = eccentricitiesForRetinaMeridianInEye(...
            radialEcc, retinaMeridian, whichEye)

    temporalRetinaMeridianName = WatsonRGCModel.enumeratedMeridianNames{1};
    nasalRetinaMeridianName = WatsonRGCModel.enumeratedMeridianNames{3};
    superiorRetinaMeridianName = WatsonRGCModel.enumeratedMeridianNames{2};
    inferiorRetinaMeridianName = WatsonRGCModel.enumeratedMeridianNames{4};
    
    % Check validity of inputs
    assert(ismember(retinaMeridian, WatsonRGCModel.enumeratedMeridianNames), ...
        sprintf('Retinal meridian name (''%s'') is invalid. It must be one of the following:\n\t''%s''\n\t''%s''\n\t''%s''\n\t''%s''', ...
        retinaMeridian, ...
        temporalRetinaMeridianName, nasalRetinaMeridianName, ...
        superiorRetinaMeridianName, inferiorRetinaMeridianName));

    assert(ismember(whichEye, {'right eye', 'left eye'}), ...
        sprintf('Eye (''%s'') is invalid. It must be either ''left eye'', or ''right eye''.', whichEye));

    assert(all(radialEcc>=0), 'radialEcc must be a non-negative scalar, vector or matrix');

    switch (retinaMeridian)
        case temporalRetinaMeridianName
            if (strcmp(whichEye, 'left eye'))
                horizontalEcc = radialEcc;
            else
                horizontalEcc = -radialEcc;
            end
            verticalEcc = abs(horizontalEcc)*0;

        case nasalRetinaMeridianName
            if (strcmp(whichEye, 'left eye'))
                horizontalEcc = -radialEcc;
            else
                horizontalEcc = radialEcc;
            end
            verticalEcc = abs(horizontalEcc)*0;

        case inferiorRetinaMeridianName
            verticalEcc = -radialEcc;
            horizontalEcc = abs(verticalEcc)*0;

        case superiodRetinaMeridianName
            verticalEcc = radialEcc;
            horizontalEcc = abs(verticalEcc)*0;
    end
end

