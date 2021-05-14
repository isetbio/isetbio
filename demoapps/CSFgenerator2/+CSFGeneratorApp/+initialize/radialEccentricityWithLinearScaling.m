function radialEccentricityWithLinearScaling(app)
    radialEccentricitySliderTicksNum = 10;
    radialEccentricitySliderMin = 0.5;
    radialEccentricitySliderMax = 40;
    radialEccSliderTickValues = round(logspace(log10(radialEccentricitySliderMin),log10(radialEccentricitySliderMax), radialEccentricitySliderTicksNum)*100)/100;
    radialEccSliderTickValues = round(linspace(0, radialEccSliderTickValues(end-1),radialEccentricitySliderTicksNum));
    
    app.roiRadialEccentricitySlider.Limits = [0 radialEccentricitySliderTicksNum-1];
    app.roiRadialEccentricitySlider.MajorTicks = 0:(radialEccentricitySliderTicksNum-1);
    app.roiRadialEccentricitySlider.MinorTicks = [];
    for k = 1:numel(app.roiRadialEccentricitySlider.MajorTicks)
        app.roiRadialEccentricitySlider.MajorTickLabels{k} = sprintf('%2.0f', radialEccSliderTickValues(k));
    end

    app.roiParams.maxEcc = radialEccSliderTickValues(end);
    CSFGeneratorApp.decode.roiRadialEccentricitySlider(app, 'valueToSlider', app.roiParams.radialEccentricityDegs);
end