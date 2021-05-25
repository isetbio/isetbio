function stimulusMosaicCenteredCheckBox(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.stimulusMosaicCenteredCheckBox.Value = value;
        case 'sliderToValue'
            app.stimParams.mosaicCenteredPosition = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
    
    % Set the 'opticalImagePositionDegs' property of the cMosaic
    if (value) 
        app.components.coneMosaic.opticalImagePositionDegs = 'mosaic-centered';
    else
        app.components.coneMosaic.opticalImagePositionDegs = app.stimParams.positionDegs;
    end
end