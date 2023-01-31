function generateR2VFTobjects(mosaicCenterParams, rfModelParams, opticsParams, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('RTVobjIndicesToBeComputed', 'all',  @(x)(((ischar(x))&&(strcmp(x, 'all'))) || isnumeric(x)));
    p.addParameter('updateRTVFobjectAtPosition', [], @(x)(isempty(x) || (numel(x)==2)));
    p.addParameter('updateRTVFobjectWithCenterConesNum', [], @(x)(isempty(x) || (numel(x)>=1)));
    p.addParameter('updateRTVFobjectWithCenterConeType', [], @(x)(isempty(x) || ismember(x, {'L', 'M'})));
    p.parse(varargin{:});

    RTVobjIndicesToBeComputed = p.Results.RTVobjIndicesToBeComputed;
    updateRTVFobjectAtPosition = p.Results.updateRTVFobjectAtPosition;
    updateRTVFobjectWithCenterConesNum = p.Results.updateRTVFobjectWithCenterConesNum;
    updateRTVFobjectWithCenterConeType = p.Results.updateRTVFobjectWithCenterConeType;

    % Generate mosaic filename and directory
    [mosaicFileName, mosaicDirectory] = midgetRGCMosaicInspector.mosaicFileName(mosaicCenterParams);

    fprintf('Loading the midget RGC mosaic from %s.\nPlease wait ...\n', mosaicFileName);
    % Load the center-connected mosaic
    load(mosaicFileName, 'theMidgetRGCmosaic');
    
    % Generate the eccentricitySamplingGrid
    eccentricitySamplingGrid = midgetRGCMosaic.eccentricitySamplingGridCoords(...
        mosaicCenterParams.positionDegs, mosaicCenterParams.sizeDegs, ...
        rfModelParams.eccentricitySamplingGridHalfSamplesNum, ...
        'hexagonal', true);

    if (~isempty(updateRTVFobjectAtPosition))
        targetPosition = [updateRTVFobjectAtPosition(1) updateRTVFobjectAtPosition(2)];
        d = sum((bsxfun(@minus, eccentricitySamplingGrid, targetPosition)).^2,2);
        [~,idx] = min(d);
        eccentricitySamplingGrid = eccentricitySamplingGrid(idx,:);
        fprintf(2, 'Will update object(s) of the RTVFobjList near (%2.2f,%2.2f) degs', ...
            eccentricitySamplingGrid(1), eccentricitySamplingGrid(2));
    end

    fprintf('\nVisualizing the mosaic. Please wait ...')
    % Visualize the mosaic
    theMidgetRGCmosaic.visualize( ...
        'eccentricitySamplingGrid', eccentricitySamplingGrid, ...
        'inputPoolingVisualization', 'centerOnly');
    fprintf('\nDone ! \n');

    % Assemble R2CVFT filename
    R2VFTobjFileName = midgetRGCMosaicInspector.R2VFTobjFileName(mosaicFileName, opticsParams, rfModelParams.H1cellIndex);

    % multiStartsNum: select from:
    % - 1 (Single start run, fastest results), 
    % - some number (Multi-start), or 
    % - inf (Global search)
    fitParams = struct();
    fitParams.multiStartsNumDoGFit = 128;

    % Where to save the fitted RTVFobjects
    fitParams.exportsDirectory = mosaicDirectory;

    tStart = cputime;

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
    
    

    if (isempty(updateRTVFobjectAtPosition))
        
        % Update fitParams with initialRetinalConePoolingParamsStruct
        fitParams.initialRetinalConePoolingParamsStruct = initialRetinalConePoolingParamsStruct;
    
        % Generate list of RTVT objects
        [theRTFVTobjList, theOpticsPositionGrid, ...
         theConesNumPooledByTheRFcenterGrid, ...
         theVisualSTFSurroundToCenterRcRatioGrid, ...
         theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid] = midgetRGCMosaic.R2VFTobjects(...
                    RTVobjIndicesToBeComputed, ...
                    theMidgetRGCmosaic, eccentricitySamplingGrid, ...
                    rfModelParams, opticsParams, fitParams);
    else
    
        fitParams.initialRetinalConePoolingParamsStruct = initialRetinalConePoolingParamsStruct;
        
        % Generate list of updated RTVT objects
        [theUpdatedRTFVTobjList, theUpdatedOpticsPositionGrid, ...
         theUpdatedConesNumPooledByTheRFcenterGrid] = midgetRGCMosaic.R2VFTobjects(...
                    RTVobjIndicesToBeComputed, ...
                    theMidgetRGCmosaic, eccentricitySamplingGrid, ...
                    rfModelParams, opticsParams, fitParams);
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

            % Query user whether to update the list of RVFTobj
            if (isempty(updateRTVFobjectWithCenterConeType))
                prompt = sprintf('Update the RTVFlist{%d}  (%d-cone RF center, BOTH L- & M-cone models - you will be asked separately for L- and M-cone model) at position (degs): (%2.2f %2.2f))? [y/n] : ', ...
                    destinationRTVFobjectIndex, theConesNumPooledByTheRFcenterGrid(destinationRTVFobjectIndex), ...
                    theOpticsPositionGrid(destinationRTVFobjectIndex,1), theOpticsPositionGrid(destinationRTVFobjectIndex,2));
            elseif strcmpi(updateRTVFobjectWithCenterConeType, 'l')
                prompt = sprintf('Update the RTVFlist{%d}  (%d-cone RF center, ONLY the L-cone model) at position (degs): (%2.2f %2.2f))? [y/n] : ', ...
                    destinationRTVFobjectIndex,  theConesNumPooledByTheRFcenterGrid(destinationRTVFobjectIndex), ...
                    theOpticsPositionGrid(destinationRTVFobjectIndex,1), theOpticsPositionGrid(destinationRTVFobjectIndex,2));
            elseif strcmpi(updateRTVFobjectWithCenterConeType, 'm')
                prompt = sprintf('Update the RTVFlist{%d}  (%d-cone RF center, ONLY the M-cone model) at position (degs): (%2.2f %2.2f))? [y/n] : ', ...
                    destinationRTVFobjectIndex, theConesNumPooledByTheRFcenterGrid(destinationRTVFobjectIndex), ...
                    theOpticsPositionGrid(destinationRTVFobjectIndex,1), theOpticsPositionGrid(destinationRTVFobjectIndex,2));
            end

            txt = lower(input(prompt,'s'));
            if isempty(txt)
                txt = 'n';
            end
            if (strcmp(txt, 'n'))
                fprintf('Will skip overwriting previous data.\n');
                continue;
            end

            reWriteRTVFTfile = true;

            % Update !
            if (isempty(updateRTVFobjectWithCenterConeType))

                % Both L- and M-. 
                
                % Ask the user about the L-cone model.
                prompt = sprintf('Update the L-cone model ? [y/n] : ');
                txt = lower(input(prompt,'s'));
                if isempty(txt)
                   txt = 'n';
                end
                if (strcmp(txt, 'y'))   
                    fprintf('Updating RTVFTobj data for the L-cone model.\n');
                    theRTFVTobjList{destinationRTVFobjectIndex}.overwriteLconeRFcomputeStruct(...
                        theUpdatedRTFVTobjList{sourceRTVFobjectIndex}.LconeRFcomputeStruct);
                end

                % Ask the user about the M-cone model.
                prompt = sprintf('Update the M-cone model ? [y/n] : ');
                txt = lower(input(prompt,'s'));
                if isempty(txt)
                   txt = 'n';
                end
                if (strcmp(txt, 'y'))
                    fprintf('Updating RTVFTobj data for the M-cone model.\n');
                    theRTFVTobjList{destinationRTVFobjectIndex}.overwriteMconeRFcomputeStruct(...
                        theUpdatedRTFVTobjList{sourceRTVFobjectIndex}.MconeRFcomputeStruct);
                end


            elseif strcmpi(updateRTVFobjectWithCenterConeType, 'l')
                fprintf('Updating RTVFTobj data for the L-cone model.\n');
                theRTFVTobjList{destinationRTVFobjectIndex}.overwriteLconeRFcomputeStruct(...
                        theUpdatedRTFVTobjList{sourceRTVFobjectIndex}.LconeRFcomputeStruct);

            elseif strcmpi(updateRTVFobjectWithCenterConeType, 'm')
                fprintf('Updating RTVFTobj data for the M-cone model.\n');
                theRTFVTobjList{destinationRTVFobjectIndex}.overwriteMconeRFcomputeStruct(...
                        theUpdatedRTFVTobjList{sourceRTVFobjectIndex}.MconeRFcomputeStruct);
            end

        end

        if (reWriteRTVFTfile)
            % Save the updated RTVFT list
            save(R2VFTobjFileName, ...
                'theRTFVTobjList', ...
                'theOpticsPositionGrid', 'theConesNumPooledByTheRFcenterGrid', ...
                'theVisualSTFSurroundToCenterRcRatioGrid', ...
                'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid', ...
                '-v7.3');
        end

    else
        % Save the computed RTVFT list
        save(R2VFTobjFileName, 'theRTFVTobjList', 'theOpticsPositionGrid', ...
                    'theConesNumPooledByTheRFcenterGrid', ...
                    'theVisualSTFSurroundToCenterRcRatioGrid', ...
                    'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid', ...
                    '-v7.3');
        fprintf('Computed R2VFTobjects saved in: %s\n', R2VFTobjFileName);
    end
end