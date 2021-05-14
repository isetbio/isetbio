function coneMosaicEccVaryingConeApertureCheckBox(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.coneMosaicEccVaryingConeApertureCheckBox.Value = value;
        case 'sliderToValue'
            app.coneMosaicParams.eccVaryingConeAperture = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end