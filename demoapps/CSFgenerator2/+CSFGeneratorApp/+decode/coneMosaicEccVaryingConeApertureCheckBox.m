function coneMosaicEccVaryingConeApertureCheckBox(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.coneMosaicEccVaryingConeApertureCheckBox.Value = value;
        case 'sliderToValue'
            app.coneMosaicParams.eccVaryingConeAperture = value;
            % Also update the current coneMosaic object's property in case
            % we are not generating a new mosaic
            app.components.coneMosaic.eccVaryingConeAperture = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end