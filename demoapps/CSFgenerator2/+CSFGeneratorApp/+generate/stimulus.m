function stimulus(app, dialog)

    deleteProgressBar = isempty(dialog);
    if (deleteProgressBar)
        % Open progressbar
        dialogBox = uiprogressdlg(app.mainView,'Title','Please Wait',...
            'Message','Generating stimulus ...');
        dialogBox.Value = 0.2;
    end

    % Generate stimParams struct fot gratingSceneEngine
    sParams = CSFGeneratorApp.generate.stimParamsStructForGratingSceneEngine(app, app.stimParams.spatialFrequencyCPD);

    % Generate the stimulus scene engine
    [app.products.demoStimulusSceneSequence, app.products.nullStimulusScene, statusReport] = ...
        CSFGeneratorApp.generate.gratingSceneEngine(sParams, []);

    % Visualize the stimulus
    CSFGeneratorApp.render.stimulusView(app, 'update');

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
        statusReportText = sprintf('Stimulus is realizable. Stimulus sequence contains %d frame(s).', numel(app.products.demoStimulusSceneSequence));
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
