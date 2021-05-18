function theOptics = optics(app, dialog)

    if (isa(app, 'ISETBioCSFGenerator'))
        appCall = true;
    else
        appCall = false;
    end
    
    if (appCall)
        deleteProgressBar = isempty(dialog);
        if (deleteProgressBar)
            % Open progressbar
            dialogBox = uiprogressdlg(app.mainView,'Title','Please Wait',...
                        'Message','Generating optics ...');
            dialogBox.Value = 0.2; 
        end
    end
    
    switch (app.opticsParams.subjectDataset)
        case 'Polans2015'
            zernikeDataBase = 'Polans2015';
        case 'Artal'
            zernikeDataBase = 'Artal';
        otherwise
            error('Unknown optics dataset: ''%s''.', app.opticsParams.subjectDataset);
    end
    
    % Generate optics appropriate for the mosaic's eccentricity
    [oiEnsemble, psfEnsemble] = ...
            app.components.coneMosaic.oiEnsembleGenerate(app.coneMosaicParams.eccentricityDegs, ...
            'zernikeDataBase', zernikeDataBase, ...
            'subjectID', app.opticsParams.subjectID, ...
            'pupilDiameterMM', app.opticsParams.pupilDiameterMM, ...
            'subtractCentralRefraction', app.opticsParams.subtractCentralRefraction);
           
     % Update all components with new optics
     app.components.psf = psfEnsemble{1};
     app.components.optics = oiEnsemble{1};
            
     if (~appCall)
        theOptics = app.components.optics;
        return;
     end
    
     % Visualize the optics
     CSFGeneratorApp.render.opticsView(app, 'update');
     
     % Update status for 'optics'
    app.statusMessages('optics')  = struct(...
        'text', sprintf('Generated optics.'), ...
        'fontColor', app.colors('good message foreground'), ...
        'backgroundColor', app.colors('good message background'), ...
        'fontWeight', 'normal');
    
    % Render the status on the status field of tab B
    CSFGeneratorApp.render.statusField(app,'B', 'optics'); 
    
    if (deleteProgressBar)
        % Close progressbar
        close(dialogBox);
    end
end
