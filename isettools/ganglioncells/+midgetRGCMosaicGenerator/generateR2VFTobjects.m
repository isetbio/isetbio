function generateR2VFTobjects(mosaicCenterParams, mosaicSurroundParams, opticsParams, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('updateRTVFobjectAtPosition', [], @(x)(isempty(x) || (numel(x)==2)));
    p.addParameter('updateRTVFobjectWithCenterConesNum', [], @(x)(isempty(x) || (numel(x)>=1)));
    p.parse(varargin{:});

    updateRTVFobjectAtPosition = p.Results.updateRTVFobjectAtPosition;
    updateRTVFobjectWithCenterConesNum = p.Results.updateRTVFobjectWithCenterConesNum;

    % Generate mosaic filename and directory
    [mosaicFileName, mosaicDirectory] = midgetRGCMosaicInspector.mosaicFileName(mosaicCenterParams);

    fprintf('Loading the midget RGC mosaic from %s.\nPlease wait ...\n', mosaicFileName);
    % Load the center-connected mosaic
    load(mosaicFileName, 'theMidgetRGCmosaic');
    
    % Generate the eccentricitySamplingGrid
    eccentricitySamplingGrid = midgetRGCMosaic.eccentricitySamplingGridCoords(...
        mosaicCenterParams.positionDegs, mosaicCenterParams.sizeDegs, ...
        mosaicSurroundParams.eccentricitySamplingGridHalfSamplesNum, ...
        'hexagonal', true);

    if (~isempty(updateRTVFobjectAtPosition))
        targetPosition = [updateRTVFobjectAtPosition(1) updateRTVFobjectAtPosition(2)];
        d = sum((bsxfun(@minus, eccentricitySamplingGrid, targetPosition)).^2,2);
        [~,idx] = min(d);
        eccentricitySamplingGrid = eccentricitySamplingGrid(idx,:);
        fprintf(2, 'Will update object(s) of the RTVFobjList near pos (degs) = (%2.2f,%2.2f)', ...
            eccentricitySamplingGrid(1), eccentricitySamplingGrid(2));
    end

    fprintf('\nVisualizing the mosaic. Please wait ...')
    % Visualize the mosaic
    theMidgetRGCmosaic.visualize( ...
        'eccentricitySamplingGrid', eccentricitySamplingGrid, ...
        'inputPoolingVisualization', 'centerOnly');
    fprintf('\nDone ! \n');

    % Assemble R2CVFT filename
    R2VFTobjFileName = midgetRGCMosaicInspector.R2VFTobjFileName(mosaicFileName, mosaicSurroundParams.H1cellIndex);

    % multiStartsNum: select from:
    % - 1 (Single start run, fastest results), 
    % - some number (Multi-start), or 
    % - inf (Global search)
    fitParams = struct();
    fitParams.multiStartsNumDoGFit = 128;

    % Where to save the fitted RTVFobjects
    fitParams.exportsDirectory = mosaicDirectory;

    tStart = cputime;

    if (isempty(updateRTVFobjectAtPosition))
        % Ask the user if he wants to use a dictionary with previously
        % fitted params to use as initial values
        usePreviousFittedParamsValues = input('Use initial values from a previous fit ? [y = YES] ', 's');
        initialRetinalConePoolingParamsStruct = [];
        if strcmpi(usePreviousFittedParamsValues, 'y')
            dropboxDir = midgetRGCMosaicInspector.localDropboxPath();
            [file,path] = uigetfile(fullfile(dropboxDir, '*.mat'), ...
                                'Select a file');
        
            if (file ~= 0)
                fName = fullfile(path,file);
                load(fName, 'retinalConePoolingParamsDictionary', 'theConesNumPooledByTheRFcenterGrid', 'theOpticsPositionGrid');
                initialRetinalConePoolingParamsStruct.dictionary = retinalConePoolingParamsDictionary;
                initialRetinalConePoolingParamsStruct.eccentricitySamplingGrid = theOpticsPositionGrid;
                initialRetinalConePoolingParamsStruct.centerConesNumGrid = theConesNumPooledByTheRFcenterGrid;
            end
        end
        
        % Update fitParams with initialRetinalConePoolingParamsStruct
        fitParams.initialRetinalConePoolingParamsStruct = initialRetinalConePoolingParamsStruct;

        % Generate list of RTVT objects
        [theRTFVTobjList, theOpticsPositionGrid, ...
         theConesNumPooledByTheRFcenterGrid, ...
         theVisualSTFSurroundToCenterRcRatioGrid, ...
         theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid] = midgetRGCMosaic.R2VFTobjects(...
                    theMidgetRGCmosaic, eccentricitySamplingGrid, ...
                    mosaicSurroundParams, opticsParams, fitParams);
    else
        fitParams.initialRetinalConePoolingParamsStruct  = [];
        
        % Generate list of updated RTVT objects
        [theUpdatedRTFVTobjList, theUpdatedOpticsPositionGrid, ...
         theUpdatedConesNumPooledByTheRFcenterGrid] = midgetRGCMosaic.R2VFTobjects(...
                    theMidgetRGCmosaic, eccentricitySamplingGrid, ...
                    mosaicSurroundParams, opticsParams, fitParams);
    end

    timeLapsedMinutes = (cputime - tStart)/60;
    fprintf('\n\n midgetRGCMosaic.R2VFTobjects were generated in %d positions and fitting took %f minutes\n', ...
            size(eccentricitySamplingGrid,1), timeLapsedMinutes);
        
    if (~isempty(updateRTVFobjectAtPosition))
        % Load previously generated theRTFVTobjList
        load(R2VFTobjFileName, 'theRTFVTobjList', ...
            'theOpticsPositionGrid', 'theConesNumPooledByTheRFcenterGrid', ...
            'theVisualSTFSurroundToCenterRcRatioGrid', ...
            'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid');

        reWriteRTVFTfile = false;
        for iUpdatedIndex = 1:numel(updateRTVFobjectWithCenterConesNum)

            % Compute sourceRTVFobjectIndex
            sourceRTVFobjectIndex = find(theUpdatedConesNumPooledByTheRFcenterGrid == updateRTVFobjectWithCenterConesNum(iUpdatedIndex));

            % Compute destination RTVFobject index
            destinationRTVFobjectIndex = midgetRGCMosaicGenerator.indexOfRTVFobject(...
                updateRTVFobjectWithCenterConesNum(iUpdatedIndex), ...
                targetPosition, ...
                theConesNumPooledByTheRFcenterGrid, ...
                theOpticsPositionGrid);

            % Query user whether to updaste the list of RVFTobj
            prompt = sprintf('Update the RTVFlist{%d}  (%d center cones at position (degs): (%2.2f %2.2f))? [y/n] : ', ...
                destinationRTVFobjectIndex, theConesNumPooledByTheRFcenterGrid(destinationRTVFobjectIndex), ...
                theOpticsPositionGrid(destinationRTVFobjectIndex,1), theOpticsPositionGrid(destinationRTVFobjectIndex,2));
            txt = lower(input(prompt,'s'));
            if isempty(txt)
                txt = 'n';
            end
            if (strcmp(txt, 'n'))
                fprintf('Will skip overwriting previous data.\n');
                continue;
            end

            reWriteRTVFTfile = true;
            fprintf('Updating RTVFTobj data.\n');
            
            % Update !
            theRTFVTobjList{destinationRTVFobjectIndex} = theUpdatedRTFVTobjList{sourceRTVFobjectIndex};
        end

        if (reWriteRTVFTfile)
            % Save the updated RTVFT list
            save(R2VFTobjFileName, ...
                'theRTFVTobjList', ...
                'theOpticsPositionGrid', 'theConesNumPooledByTheRFcenterGrid', ...
                'theVisualSTFSurroundToCenterRcRatioGrid', ...
                'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid');
        end

    else
        % Save the computed RTVFT list
        save(R2VFTobjFileName, 'theRTFVTobjList', 'theOpticsPositionGrid', ...
                    'theConesNumPooledByTheRFcenterGrid', ...
                    'theVisualSTFSurroundToCenterRcRatioGrid', ...
                    'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid');
        fprintf('Computed R2VFTobjects saved in: %s\n', R2VFTobjFileName);
    end
end