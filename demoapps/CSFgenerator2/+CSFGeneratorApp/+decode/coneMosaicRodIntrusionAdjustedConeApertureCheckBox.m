function coneMosaicRodIntrusionAdjustedConeApertureCheckBox(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.coneMosaicRodIntrusionAdjustedConeApertureCheckBox.Value = value;
        case 'sliderToValue'
            app.coneMosaicParams.rodIntrusionAdjustedConeAperture = value;
            % Also update the current coneMosaic object's property in case
            % we are not generating a new mosaic
            app.components.coneMosaic.rodIntrusionAdjustedConeAperture = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end