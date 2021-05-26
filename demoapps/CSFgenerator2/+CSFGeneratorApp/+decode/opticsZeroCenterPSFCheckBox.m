function opticsZeroCenterPSFCheckBox(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.opticsZeroCenterPSFCheckBox.Value = value;
        case 'sliderToValue'
            app.opticsParams.zeroCenterPSF = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end