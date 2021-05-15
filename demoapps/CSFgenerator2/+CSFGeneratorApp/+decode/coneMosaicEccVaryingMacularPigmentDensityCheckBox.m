function coneMosaicEccVaryingMacularPigmentDensityCheckBox(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.coneMosaicEccVaryingMacularPigmentDensityCheckBox.Value = value;
        case 'sliderToValue'
            app.coneMosaicParams.eccVaryingMacularPigmentDensity = value;
            % Also update the current coneMosaic object's property in case
            % we are not generating a new mosaic
            app.components.coneMosaic.eccVaryingMacularPigmentDensity = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end