function stimulusSpatialPositionXSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.stimulusSpatialPositionXSpinner.Value = value;
        case 'sliderToValue'
            app.stimParams.positionDegs(1) = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end