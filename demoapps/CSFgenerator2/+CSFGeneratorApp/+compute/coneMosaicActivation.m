function coneMosaicActivation(app, dialog)

    deleteProgressBar = isempty(dialog);
    if (deleteProgressBar)
        % Open progressbar
        dialogBox = uiprogressdlg(app.mainView,'Title','Please Wait',...
                    'Message','Computing cone mosaic activation ...');
        dialogBox.Value = 0.2; 
    end
    
    
    % Check where to position the optical image
    if (app.stimParams.mosaicCenteredPosition)
        opticalImagePositionDegs = 'mosaic-centered';
    else
        opticalImagePositionDegs = app.stimParams.positionDegs;
    end
       
    % Compute the activation to the null stimulus 
    nullActivation = app.components.coneMosaic.compute(...
        oiCompute(app.products.nullStimulusScene, app.components.optics));
    
    % Compute response to first stimulus frame only
    theScene = app.products.demoStimulusSceneSequence{1};
    
    % Compute nInstances noisy response instances
    nInstances = 4;
    [app.products.noiseFreeConeMosaicActivation, ...
     app.products.noisyConeMosaicActivationInstances] = app.components.coneMosaic.compute(...
        oiCompute(theScene, app.components.optics), ...
        'nTrials', nInstances);
    
    % Compute noise-free modulation 
    app.products.noiseFreeConeMosaicModulation = ...
        100*(app.products.noiseFreeConeMosaicActivation - nullActivation)./nullActivation;
    
    % Compute noisy modulations
    app.products.noisyConeMosaicModulationInstances = ...
        100*(app.products.noisyConeMosaicActivationInstances - nullActivation)./nullActivation;
   
    % Compute the residual activation
    if (isfield(app.products,'lastNoiseFreeConeMosaicActivation')) && ...
       (numel(app.products.lastNoiseFreeConeMosaicActivation) == numel(app.products.noiseFreeConeMosaicActivation))
        app.products.residualConeMosaicActivation = 100 * ...
            (app.products.noiseFreeConeMosaicActivation-app.products.lastNoiseFreeConeMosaicActivation) ./ ...
             app.products.noiseFreeConeMosaicActivation;
    else
        app.products.residualConeMosaicActivation = 0*app.products.noiseFreeConeMosaicActivation;
    end
    
    % Save the last mean activation
    app.products.lastNoiseFreeConeMosaicActivation = app.products.noiseFreeConeMosaicActivation;


    if (deleteProgressBar)
        % Close progressbar
        close(dialogBox);
    end
    
end
