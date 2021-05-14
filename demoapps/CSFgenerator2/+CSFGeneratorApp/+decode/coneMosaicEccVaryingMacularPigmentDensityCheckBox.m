function coneMosaicEccVaryingMacularPigmentDensityCheckBox(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.coneMosaicEccVaryingMacularPigmentDensityCheckBox.Value = value;
        case 'sliderToValue'
            app.coneMosaicParams.eccVaryingMacularPigmentDensity = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end