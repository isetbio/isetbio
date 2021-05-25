function stimulusSpatialPositionYSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.stimulusSpatialPositionYSpinner.Value = value;
        case 'sliderToValue'
            app.stimParams.positionDegs(2) = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
    
    % Update the y-position of the 'opticalImagePositionDegs' property of the cMosaic
    if (~app.stimParams.mosaicCenteredPosition) 
        app.components.coneMosaic.opticalImagePositionDegs(2) = value;
    end
end