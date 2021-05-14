function coneMosaicActivationTypeSwitch(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.coneMosaicActivationTypeSwitch.Value = value;
        case 'sliderToValue'
            app.viewModes.coneMosaicActivationType = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end
