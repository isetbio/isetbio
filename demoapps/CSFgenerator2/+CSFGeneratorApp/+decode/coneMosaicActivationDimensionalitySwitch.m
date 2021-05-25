function coneMosaicActivationDimensionalitySwitch(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.coneMosaicActivationDimensionalitySwitch.Value = value;
        case 'sliderToValue'
            app.viewModes.coneMosaicActivationDimensionality = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end