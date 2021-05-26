function opticsFlipPSFUpsideDownCheckBox(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.opticsFlipPSFUpsideDownCheckBox.Value = value;
        case 'sliderToValue'
            app.opticsParams.flipPSFUpsideDown = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end