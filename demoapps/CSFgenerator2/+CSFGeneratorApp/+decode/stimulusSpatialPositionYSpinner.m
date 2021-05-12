function stimulusSpatialPositionYSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.stimulusSpatialPositionYSpinner.Value = value;
        case 'sliderToValue'
            app.stimParams.positionDegs(2) = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end