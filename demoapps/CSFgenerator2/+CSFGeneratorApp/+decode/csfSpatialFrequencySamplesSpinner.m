function csfSpatialFrequencySamplesSpinner(app, direction, value)
    switch direction
        case 'valueToSlider'
            app.csfSpatialFrequencySamplesSpinner.Value = value;
        case 'sliderToValue'
           app.csfParams.spatialFrequencySamples = value;
           CSFGeneratorApp.render.csfView(app, 'update');
        otherwise
            error('Unknown direction; ''%s''.'\n', direction);
    end
    
    % Colormap for the different spatial frequencies
    app.csfLineColors = brewermap(value+2, 'Spectral');

end