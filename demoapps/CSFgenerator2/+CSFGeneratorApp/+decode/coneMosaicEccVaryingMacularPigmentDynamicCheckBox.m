function coneMosaicEccVaryingMacularPigmentDynamicCheckBox(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.coneMosaicEccVaryingMacularPigmentDynamicCheckBox.Value = value;
        case 'sliderToValue'
            app.coneMosaicParams.eccVaryingMacularPigmentDynamic = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end