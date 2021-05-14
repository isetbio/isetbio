function coneMosaicActivationSignalSwitch(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.coneMosaicActivationSignalSwitch.Value = value;
        case 'sliderToValue'
            app.viewModes.coneMosaicActivationSignal = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end