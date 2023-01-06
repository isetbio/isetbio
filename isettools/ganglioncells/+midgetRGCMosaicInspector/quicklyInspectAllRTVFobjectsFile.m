function quicklyInspectAllRTVFobjectsFile()
    
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


    for iObj = 1:numel(theRTFVTobjList)
        waitbar(0.3+(iObj/numel(theRTFVTobjList))*0.7,progressBar, sprintf('Processing R2VF obj %d of %d', iObj, numel(theRTFVTobjList)));
        pause(0.1);
        midgetRGCMosaicInspector.peekIntoRTVFobj(theRTFVTobjList{iObj}, iObj, theOpticsPositionGrid, theConesNumPooledByTheRFcenterGrid, 1000+iObj);
    end

    close(progressBar);
end