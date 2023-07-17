function [horizontalRetinalMeridian, verticalRetinalMeridian] = retinalMeridiansForEccentricityInEye(horizontalEcc, verticalEcc, whichEye)

    temporalRetinaMeridianName = WatsonRGCModel.enumeratedMeridianNames{1};
    nasalRetinaMeridianName = WatsonRGCModel.enumeratedMeridianNames{3};
    superiorRetinaMeridianName = WatsonRGCModel.enumeratedMeridianNames{2};
    inferiorRetinaMeridianName = WatsonRGCModel.enumeratedMeridianNames{4};
   
    assert(ismember(whichEye, {'right eye', 'left eye'}), ...
        sprintf('Eye (''%s'') is invalid. It must be either ''left eye'', or ''right eye''.', whichEye));

    switch (whichEye)
        case 'left eye'
            if (horizontalEcc >= 0)
                horizontalRetinalMeridian = temporalRetinaMeridianName;
            else
                horizontalRetinalMeridian = nasalRetinaMeridianName;
            end

        case 'right eye'
            if (horizontalEcc >= 0)
                horizontalRetinalMeridian = nasalRetinaMeridianName;
            else
                horizontalRetinalMeridian = temporalRetinaMeridianName;
            end
    end

    if (verticalEcc >= 0)
        verticalRetinalMeridian = superiorRetinaMeridianName;
    else
        verticalRetinalMeridian = inferiorRetinaMeridianName;
    end

end