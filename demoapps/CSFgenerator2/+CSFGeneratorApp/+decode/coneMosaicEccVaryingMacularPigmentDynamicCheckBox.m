function coneMosaicEccVaryingMacularPigmentDynamicCheckBox(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.coneMosaicEccVaryingMacularPigmentDynamicCheckBox.Value = value;
        case 'sliderToValue'
            app.coneMosaicParams.eccVaryingMacularPigmentDynamic = value;
            % Also update the current coneMosaic object's property in case
            % we are not generating a new mosaic
            app.components.coneMosaic.eccVaryingMacularPigmentDensityDynamic = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end