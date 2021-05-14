function coneMosaicEccVaryingConeApertureBlurCheckBox(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.coneMosaicEccVaryingConeApertureBlurCheckBox.Value = value;
        case 'sliderToValue'
            app.coneMosaicParams.eccVaryingConeApertureBlur = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end