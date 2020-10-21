function visualMeridianName = retinalMeridianNameToVisualMeridianName(retinalMeridianName, whichEye)
% Convert retinal meridian name at a specified eye to the visual meridian name
% of the right eye, which is the meridian used by Watson

    % Validate the meridian name
    RGCmodels.Watson.validate.meridianName(retinalMeridianName);
    
    % Validate the eye name
    RGCmodels.Watson.validate.eyeName(whichEye);
    
    % Get the meridian names
    nasalMeridian = RGCmodels.Watson.constants.nasalMeridian;
    temporalMeridian = RGCmodels.Watson.constants.temporalMeridian;
    superiorMeridian = RGCmodels.Watson.constants.superiorMeridian;
    inferiorMeridian = RGCmodels.Watson.constants.inferiorMeridian;
    
    % Make dictionary mapping retinal meridian name -> visual meridian name
    visualMeridianNames = containers.Map();
    if (strcmp(whichEye, RGCmodels.Watson.constants.leftEye))
        visualMeridianNames(nasalMeridian) = nasalMeridian;
        visualMeridianNames(temporalMeridian) = temporalMeridian;
    else
        visualMeridianNames(nasalMeridian) = temporalMeridian;
        visualMeridianNames(temporalMeridian) = nasalMeridian;
    end
    visualMeridianNames(superiorMeridian) = inferiorMeridian;
    visualMeridianNames(inferiorMeridian) = superiorMeridian;
    
    % Return corresponding visual meridian name
    visualMeridianName = visualMeridianNames(retinalMeridianName);
end

