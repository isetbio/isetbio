function [horizontalRetinalMeridian, verticalRetinalMeridian] = retinalMeridiansForEccentricityInEye(horizontalEcc, verticalEcc, whichEye)

    temporalRetinaMeridianName = RGCmodels.Watson.constants.indexedMeridians{1};
    nasalRetinaMeridianName = RGCmodels.Watson.constants.indexedMeridians{3};
    superiorRetinaMeridianName = RGCmodels.Watson.constants.indexedMeridians{2};
    inferiorRetinaMeridianName = RGCmodels.Watson.constants.indexedMeridians{4};
   
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