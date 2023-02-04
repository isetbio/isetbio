function exportRetinalConePoolingParamsForAllFittedRTVFTobjects()
    
    dropboxDir = midgetRGCMosaicInspector.localDropboxPath();
    midgetRGCMosaicInspector.say('Select the RTVF-objects file to export ');
    [file,path] = uigetfile(fullfile(dropboxDir, '*.mat'), ...
                        'Select the RTVF-objects file to export');

    if (file == 0)
        return;
    end

    midgetRGCMosaicInspector.say('Loading computed R2VFT objects. Please wait');

    fName = fullfile(path,file);
    load(fName, 'theRTFVTobjList', ...
        'theConesNumPooledByTheRFcenterGrid', ...
        'theOpticsPositionGrid');

    retinalConePoolingParamsDictionary = containers.Map();
    for iRTVobjIndex = 1:numel(theRTFVTobjList)
        
        theKey = sprintf('%d center cones at (%2.3f,%2.3f)', ...
            theConesNumPooledByTheRFcenterGrid(iRTVobjIndex), ...
            theOpticsPositionGrid(iRTVobjIndex,1), ...
            theOpticsPositionGrid(iRTVobjIndex,2));

        theRTVFTobj = theRTFVTobjList{iRTVobjIndex};
        s = struct();   
        s.LconeRetinalConePoolingParams = theRTVFTobj.LconeRFcomputeStruct.retinalConePoolingParams;
        s.MconeRetinalConePoolingParams = theRTVFTobj.MconeRFcomputeStruct.retinalConePoolingParams;

        retinalConePoolingParamsDictionary(theKey) = s;
    end

   
    midgetRGCMosaicInspector.say('Loaded all computed R2VFT objects. Please specify a file name for export');

    suggestedExportFileName = strrep(fName, '.mat', '_FittedRTVTparamsOnly.mat');
    [file,path] = uiputfile(suggestedExportFileName);
    fName = fullfile(path,file);

    save(fName, 'retinalConePoolingParamsDictionary', 'theConesNumPooledByTheRFcenterGrid', ...
        'theOpticsPositionGrid');

    fprintf('Fitted RTVFTparams for all objects was exported to %s\n', fName);
end