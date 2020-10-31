function meridianConversions()
% Test meridian convertions from visual to retinal space
    nasalMeridian = RGCmodels.Watson.constants.nasalMeridian;
    temporalMeridian = RGCmodels.Watson.constants.temporalMeridian;
    superiorMeridian = RGCmodels.Watson.constants.superiorMeridian;
    inferiorMeridian = RGCmodels.Watson.constants.inferiorMeridian;
    
    % Left eye tests
    theEye = RGCmodels.Watson.constants.leftEye;
    assert(strcmp(RGCmodels.Watson.convert.retinalMeridianNameToVisualMeridianName(nasalMeridian, theEye),nasalMeridian), ...
        sprintf('Incorrect conversion of nasal retinal meridian in the %s', theEye));
    assert(strcmp(RGCmodels.Watson.convert.retinalMeridianNameToVisualMeridianName(temporalMeridian, theEye),temporalMeridian), ...
        sprintf('Incorrect conversion of temporal retinal meridian in %s', theEye));
    assert(strcmp(RGCmodels.Watson.convert.retinalMeridianNameToVisualMeridianName(superiorMeridian, theEye),inferiorMeridian), ...
        sprintf('Incorrect conversion of superior retinal meridian in %s', theEye));
    assert(strcmp(RGCmodels.Watson.convert.retinalMeridianNameToVisualMeridianName(inferiorMeridian, theEye),superiorMeridian), ...
        sprintf('Incorrect conversion of inferior retinal meridian in %s', theEye));
    
    % Right eye tests
    theEye = RGCmodels.Watson.constants.rightEye;
    assert(strcmp(RGCmodels.Watson.convert.retinalMeridianNameToVisualMeridianName(nasalMeridian, theEye),temporalMeridian), ...
        sprintf('Incorrect conversion of nasal retinal meridian in the %s', theEye));
    assert(strcmp(RGCmodels.Watson.convert.retinalMeridianNameToVisualMeridianName(temporalMeridian, theEye),nasalMeridian), ...
        sprintf('Incorrect conversion of temporal retinal meridian in %s', theEye));
    assert(strcmp(RGCmodels.Watson.convert.retinalMeridianNameToVisualMeridianName(superiorMeridian, theEye),inferiorMeridian), ...
        sprintf('Incorrect conversion of superior retinal meridian in %s', theEye));
    assert(strcmp(RGCmodels.Watson.convert.retinalMeridianNameToVisualMeridianName(inferiorMeridian, theEye),superiorMeridian), ...
        sprintf('Incorrect conversion of inferior retinal meridian in %s', theEye));

end

