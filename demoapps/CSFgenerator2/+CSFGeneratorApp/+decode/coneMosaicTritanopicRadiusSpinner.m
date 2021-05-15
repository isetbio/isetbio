function coneMosaicTritanopicRadiusSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.coneMosaicTritanopicRadiusSpinner.Value = value;
        case 'sliderToValue'
           app.coneMosaicParams.tritanopicRadiusDegs = value;
           % Also update the current coneMosaic object's property in case
           % we are not generating a new mosaic
           app.components.coneMosaic.tritanopicRadiusDegs = value;
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
end