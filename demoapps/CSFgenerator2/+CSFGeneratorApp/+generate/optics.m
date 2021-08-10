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
            rankedSubjectIDs = PolansOptics.constants.subjectRanking;
            subjectID = rankedSubjectIDs(app.opticsParams.subjectRank);
            fprintf('Will generate optics for Polans subject %d (rank: %d)\n', subjectID,app.opticsParams.subjectRank);
            % Determine if we need to subtract the subject's central refraction
            subtractCentralRefraction = PolansOptics.constants.subjectRequiresCentralRefractionCorrection(subjectID);

        case 'Artal2012'
            zernikeDataBase = 'Artal2012';
            rankedSubjectIDs = ArtalOptics.constants.subjectRanking(app.roiParams.whichEye);
            subjectID = rankedSubjectIDs(app.opticsParams.subjectRank);
            fprintf('Will generate optics for Artal subject %d (rank: %d)\n', subjectID,app.opticsParams.subjectRank);
            % Determine if we need to subtract the subject's central refraction
            subtractCentralRefraction = ArtalOptics.constants.subjectRequiresCentralRefractionCorrection(app.roiParams.whichEye,subjectID);

        otherwise
            error('Unknown optics dataset: ''%s''.', app.opticsParams.subjectDataset);
    end
    
    
    % Generate optics appropriate for the mosaic's eccentricity
    [oiEnsemble, psfEnsemble] = ...
            app.components.coneMosaic.oiEnsembleGenerate(app.coneMosaicParams.eccentricityDegs, ...
            'zernikeDataBase', zernikeDataBase, ...
            'subjectID', subjectID, ...
            'pupilDiameterMM', app.opticsParams.pupilDiameterMM, ...
            'subtractCentralRefraction', subtractCentralRefraction, ...
            'zeroCenterPSF', app.opticsParams.zeroCenterPSF, ...
            'flipPSFUpsideDown', app.opticsParams.flipPSFUpsideDown, ...
            'wavefrontSpatialSamples', app.opticsParams.wavefrontSpatialSamples ...
            );
           
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
