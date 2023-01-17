function exportRetinalConePoolingParamsForAllFittedRTVFTobjects()
    
    dropboxDir = midgetRGCMosaicInspector.localDropboxPath();
    [file,path] = uigetfile(fullfile(dropboxDir, '*.mat'), ...
                        'Select an RTVF file');

    if (file == 0)
        return;
    end

    progressBar = waitbar(0.2,'Loading all computed R2VFT objects. Please wait ...');
    pause(.1);

    fName = fullfile(path,file);
    load(fName, 'theRTFVTobjList', ...
        'theConesNumPooledByTheRFcenterGrid', ...
        'theOpticsPositionGrid');

    retinalConePoolingParamsDictionary = containers.Map();
    for iRTVobjIndex = 1:numel(theRTFVTobjList)
        theRTVFTobj = theRTFVTobjList{iRTVobjIndex};
        
        theKey = sprintf('%d center cones at (%2.3f,%2.3f)', ...
            theConesNumPooledByTheRFcenterGrid(iRTVobjIndex), ...
            theOpticsPositionGrid(iRTVobjIndex,1), ...
            theOpticsPositionGrid(iRTVobjIndex,2));
        retinalConePoolingParamsDictionary(theKey) = theRTVFTobj.rfComputeStruct.retinalConePoolingParams;
    end

    close(progressBar);

    suggestedExportFileName = strrep(fName, '.mat', '_FittedRTVTparamsOnly.mat');
    [file,path] = uiputfile(suggestedExportFileName);
    fName = fullfile(path,file);

    fprintf('Fitted RTVFTparams for all objects will be exported to %s\n', fName);
    fprintf('Hit Enter to continue ...');
    pause
    save(fName, 'retinalConePoolingParamsDictionary', 'theConesNumPooledByTheRFcenterGrid', ...
        'theOpticsPositionGrid');

end