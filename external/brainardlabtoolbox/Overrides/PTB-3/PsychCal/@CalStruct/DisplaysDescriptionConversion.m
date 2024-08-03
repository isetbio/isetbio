function displayDescription = DisplaysDescriptionConversion(obj, propertyValue)
    % In the new-style cal we store the description for all attached
    % displays, in the old-style cal we stored the description of the display
    % that was calibrated. So do the same here.
    displayDescription = struct(propertyValue(obj.describe___whichScreen));
end