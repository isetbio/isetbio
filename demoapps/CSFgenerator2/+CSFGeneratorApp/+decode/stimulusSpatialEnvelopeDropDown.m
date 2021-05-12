function stimulusSpatialEnvelopeDropDown(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.stimulusSpatialEnvelopeDropDown.Value = value;
        case 'sliderToValue'
            app.stimParams.spatialEnvelope = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end