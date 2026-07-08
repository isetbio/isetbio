function rightEyeEccVisualAngles = retinalAnglesAtSpecificEyeToAnglesInRightVisualField(eccRetinalAngles, whichEye)
    xCoords = cosd(eccRetinalAngles);
    yCoords = sind(eccRetinalAngles);
    
    % Flip y-coord
    yCoords = -yCoords;
    if (strcmp(whichEye, RGCmodels.Watson.constants.rightEye))
        xCoords = -xCoords;
    end
    
    rightEyeEccVisualAngles = atan2d(yCoords, xCoords);
    rightEyeEccVisualAngles = mod(rightEyeEccVisualAngles + 360, 360);
end

