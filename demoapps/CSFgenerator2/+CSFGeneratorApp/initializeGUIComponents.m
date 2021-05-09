function initializeGUIComponents(app)

    % Initialize the visualField GUI components
    initializePolarEccentricityGUIComponents(app);
    initializeRadialEccentricityGUIComponents(app);
    initializeFieldOfViewGUIComponents(app);
    
    app.visualFieldEyeSwitch.Value = app.visualFieldParams.eye;
    app.visualFieldMagnificationFactorEditField.Value = app.visualFieldParams.magnificationFactorMicronsPerDeg;
end

function initializeFieldOfViewGUIComponents(app)
    
    app.visualFieldFieldOfViewSlider.Limits = [0.0 10];
    app.visualFieldFieldOfViewSlider.MajorTicks = 0:1:10;
    app.visualFieldFieldOfViewSlider.MajorTickLabels = {'0.1', '1', '2', '3', '4', '5', '6','7', '8', '9', '10'};
    app.visualFieldFieldOfViewSlider.MinorTicks = [];
    
    app.visualFieldFieldOfViewSlider.Value = app.visualFieldParams.fieldOfViewDegs;
    app.visualFieldMagnificationFactorEditField.Value = app.visualFieldParams.magnificationFactorMicronsPerDeg;
    
    CSFGeneratorApp.decodeFieldOfViewSlider(app, 'valueToSlider', app.visualFieldParams.fieldOfViewDegs);
end


function initializePolarEccentricityGUIComponents(app)
    app.visualFieldPolarEccentricityKnob.Limits = [0 360];
    app.visualFieldPolarEccentricityKnob.MajorTicks = 0:45:360;
    app.visualFieldPolarEccentricityKnob.MajorTickLabels = {'360', '315', '270', '225', '180', '135',  '90',  '45', '0'};
    app.visualFieldPolarEccentricityKnob.MinorTicks = 0:15:360;
    
    CSFGeneratorApp.decodePolarEccentricitySlider(app, 'valueToSlider', app.visualFieldParams.polarEccentricityDegs);
end


function initializeRadialEccentricityGUIComponents(app)

    radialEccentricitySliderTicksNum = 10;
    radialEccentricitySliderMin = 0.5;
    radialEccentricitySliderMax = 59;
    radialEccSliderTickValues = round(logspace(log10(radialEccentricitySliderMin),log10(radialEccentricitySliderMax), radialEccentricitySliderTicksNum)*100)/100;
    radialEccSliderTickValues = [0 radialEccSliderTickValues(1:end-1)];
    
    app.visualFieldRadialEccentricitySlider.Limits = [0 radialEccentricitySliderTicksNum-1];
    app.visualFieldRadialEccentricitySlider.MajorTicks = 0:(radialEccentricitySliderTicksNum-1);
    app.visualFieldRadialEccentricitySlider.MinorTicks = [];
    for k = 1:numel(app.visualFieldRadialEccentricitySlider.MajorTicks)
        if (radialEccSliderTickValues(k) < 1.0)
            app.visualFieldRadialEccentricitySlider.MajorTickLabels{k} = sprintf('%2.1f', radialEccSliderTickValues(k));
        else
            app.visualFieldRadialEccentricitySlider.MajorTickLabels{k} = sprintf('%2.0f', radialEccSliderTickValues(k));
        end
    end

    CSFGeneratorApp.decodeRadialEccentricitySlider(app, 'valueToSlider', app.visualFieldParams.radialEccentricityDegs);
end
