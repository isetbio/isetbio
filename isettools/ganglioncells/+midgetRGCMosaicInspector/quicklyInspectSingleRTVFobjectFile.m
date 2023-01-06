function quicklyInspectSingleRTVFobjectFile()

    dropboxDir = midgetRGCMosaicInspector.localDropboxPath();
    [file,path] = uigetfile(fullfile(dropboxDir, '*.mat'), ...
                        'Select an RTVF file');

    if (file == 0)
        return;
    end

    fName = fullfile(path,file);
    load(fName, 'obj', 'iRTVobjIndex', ...
        'theConesNumPooledByTheRFcenterGrid', ...
        'theSamplingPositionGrid');

    midgetRGCMosaicInspector.peekIntoRTVFobj(obj, iRTVobjIndex, theSamplingPositionGrid, theConesNumPooledByTheRFcenterGrid, 999);
end