function stimulus(app, dialog)

    deleteProgressBar = isempty(dialog);
    if (deleteProgressBar)
        % Open progressbar
        dialogBox = uiprogressdlg(app.mainView,'Title','Please Wait',...
                    'Message','Generating stimulus ...');
        dialogBox.Value = 0.2; 
    end
    
    % Form sParams struct
    switch (app.stimParams.spatialEnvelope)
        case 'disk'
            app.stimParams.spatialEnvelopeRadiusDegs = app.stimParams.sizeDegs/2;
        case 'rect'
            app.stimParams.spatialEnvelopeRadiusDegs = app.stimParams.sizeDegs/(2*sqrt(2.0));
        case 'soft'
            app.stimParams.spatialEnvelopeRadiusDegs = app.stimParams.sizeDegs/6;  
        otherwise
            error('Unknown spatial envelope: ''%s''.', app.stimParams.spatialEnvelope);
    end
            
    sParams = struct(...
                'contrast', 1.0, ...
                'chromaDir', [app.stimParams.LconeContrast app.stimParams.MconeContrast app.stimParams.SconeContrast]/100, ...
                'meanLuminanceCdPerM2', app.stimParams.meanLuminanceCdM2, ...
                'sf', app.stimParams.spatialFrequencyCPD, ...
                'fovDegs', app.stimParams.sizeDegs, ...
                'orientation', app.stimParams.orientationDegs, ...
                'spatialEnvelope', app.stimParams.spatialEnvelope, ...
                'spatialEnvelopeRadiusDegs', app.stimParams.spatialEnvelopeRadiusDegs, ...
                'spatialPhase', app.stimParams.spatialPhaseDegs, ...
                'pixelsNum', app.stimParams.resolutionPixels, ...
                'minPixelsNumPerCycle', app.stimParams.minPixelsNumPerCycle, ...
                'spectralSupport', (app.stimParams.wavelengthSupportMin : app.stimParams.wavelengthSupportStepSize: app.stimParams.wavelengthSupportMax), ...
                'presentationMode', app.stimParams.presentationMode, ...
                'duration', app.stimParams.durationSec/1000, ...
                'warningInsteadOfErrorOnOutOfGamut', true ...
     );
            
     % Call the external CSFgenerator.computeStimulusSceneEngine method to
     % generate the stimulus scene
     [app.products.demoStimulusSceneSequence, app.products.nullStimulusScene, statusReport] = ...
         CSFGeneratorApp.generate.stimulusSceneEngine(sParams, []);
            
     % Visualize the stimulus
     CSFGeneratorApp.render.stimulusView(app, 'update');
     
     dialogBox.Value = 0.8; 
     dialogBox.Message = 'Computing cone mosaic response. Please wait ...';
    
     % Compute the cone mosaic activation
     CSFGeneratorApp.compute.coneMosaicActivation(app, dialogBox);
     
     % Visualize the cone mosaic
     CSFGeneratorApp.render.coneMosaicView(app, 'update');
    
     % Update the status for the tab
     setStatusMessage(app, statusReport);
     
     if (deleteProgressBar)
        % Close progressbar
        close(dialogBox);
     end
end

function setStatusMessage(app, statusReport)

    flagNames = fieldnames(statusReport);
    if (~isempty(flagNames))
        statusReportText = '';
        for k = 1:numel(flagNames)
            statusReportText = sprintf('%s - %s', statusReportText,flagNames{k});
        end
        statusReportText = sprintf('%s -',statusReportText);
        statusFontWeight = 'Bold';
        statusFontColor = app.colors('problem message foreground');
        statusBackgroundColor = app.colors('problem message background');
    else
        statusReportText = sprintf('Stimulus is realizable. Generated stimulus sequence with %d frame(s).', numel(app.products.demoStimulusSceneSequence));
        statusFontWeight = 'Normal';
        statusFontColor = app.colors('good message foreground');
        statusBackgroundColor = app.colors('good message background');
    end
      
    % Update status for 'stimulus'
    app.statusMessages('stimulus')  = struct(...
        'text', statusReportText, ...
        'fontColor', statusFontColor, ...
        'backgroundColor', statusBackgroundColor, ...
        'fontWeight', statusFontWeight);
    
    % Render the status on the  status field of tab A
    CSFGeneratorApp.render.statusField(app,'A', 'stimulus');            
end
