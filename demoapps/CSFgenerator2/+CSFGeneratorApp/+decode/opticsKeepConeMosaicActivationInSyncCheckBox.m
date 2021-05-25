function opticsKeepConeMosaicActivationInSyncCheckBox(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.opticsKeepConeMosaicActivationInSyncCheckBox.Value = value;
        case 'sliderToValue'
            app.viewModes.opticsKeepConeMosaicActivationInSync = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end