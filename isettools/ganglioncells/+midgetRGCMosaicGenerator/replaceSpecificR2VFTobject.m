function replaceSpecificR2VFTobject()
    
    dropboxDir = midgetRGCMosaicInspector.localDropboxPath();

    fprintf(2,'****** Select the file containing the updated R2VFT object *****\n');
    % Get the replacement R2VFT file 
    [file,path] = uigetfile(fullfile(dropboxDir, '*.mat'), ...
                        'Select the file containing the updated R2VFT object');

    if (file == 0)
        return;
    end

    fName = fullfile(path,file);
    load(fName, 'obj', 'iRTVobjIndex', ...
        'theConesNumPooledByTheRFcenterGrid', ...
        'theSamplingPositionGrid');

    theTargetConesNum = theConesNumPooledByTheRFcenterGrid(iRTVobjIndex);
    theTargetSamplingPosition = theSamplingPositionGrid(iRTVobjIndex,:);
    theReplacingRTVFobj = obj;

    fprintf('The file contains the RTVFobject for %d cones at position (%2.2f, %2.2f)\n', ...
        theConesNumPooledByTheRFcenterGrid(iRTVobjIndex), ...
        theSamplingPositionGrid(iRTVobjIndex,1), theSamplingPositionGrid(iRTVobjIndex,2));

    % Get the file containing the list of all RTVF objects for this mosaic
    fprintf(2,'***** Select the file containing the list of all RTVF objects for this mosaic ****\n');
    [file,path] = uigetfile(fullfile(dropboxDir, '*.mat'), ...
                        'Select the file containing the list of all RTVF objects for this mosaic');

    if (file == 0)
        return;
    end

    progressBar = waitbar(0.2,'Loading all computed R2VFT objects. Please wait ...');
    pause(.1);

    fName = fullfile(path,file);
    load(fName, 'theRTFVTobjList', ...
        'theConesNumPooledByTheRFcenterGrid', ...
        'theOpticsPositionGrid');

    theIndexOfTheRTVFobjToBeReplaced = midgetRGCMosaicGenerator.indexOfRTVFobject(...
        theTargetConesNum, theTargetSamplingPosition, ...
        theConesNumPooledByTheRFcenterGrid, ...
        theOpticsPositionGrid);

    fprintf('The RTVF list at index %d will be replaced. If you are sure, hit ENTER to continue. \n', theIndexOfTheRTVFobjToBeReplaced);
    pause

    % Update list with replacing object
    theRTFVTobjList{theIndexOfTheRTVFobjToBeReplaced} = theReplacingRTVFobj;

    % Save new list
    save(fName, 'theRTFVTobjList', ...
        'theConesNumPooledByTheRFcenterGrid', ...
        'theOpticsPositionGrid');

    fprintf('The file containing the RTVFlist (%s) was updated\n', fName);

end