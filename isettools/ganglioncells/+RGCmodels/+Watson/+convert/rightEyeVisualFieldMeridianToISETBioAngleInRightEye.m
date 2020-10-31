function isetbioAngle = rightEyeVisualFieldMeridianToISETBioAngleInRightEye(meridianName)
% Convert the rightEye visual field meridianName to the corresponding
% ISETBio angle in the right eye

    switch (meridianName)
        case RGCmodels.Watson.constants.temporalMeridian
            isetbioAngle = 180;
        case RGCmodels.Watson.constants.superiorMeridian
            isetbioAngle = 90;
        case RGCmodels.Watson.constants.nasalMeridian
            isetbioAngle = 0;
        case RGCmodels.Watson.constants.inferiorMeridian
            isetbioAngle = 270;
    end

end

