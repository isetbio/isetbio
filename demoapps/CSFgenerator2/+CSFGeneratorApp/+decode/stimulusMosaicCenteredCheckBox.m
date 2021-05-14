function stimulusMosaicCenteredCheckBox(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.stimulusMosaicCenteredCheckBox.Value = value;
        case 'sliderToValue'
            app.stimParams.mosaicCenteredPosition = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end