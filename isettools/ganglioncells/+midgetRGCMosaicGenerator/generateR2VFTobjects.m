function generateR2VFTobjects(mosaicCenterParams, mosaicSurroundParams, opticsParams, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('updateRTVFobjectAtPosition', [], @(x)(isempty(x) || (numel(x)==2)));
    p.addParameter('updateRTVFobjectWithCenterConesNum', [], @(x)(isempty(x) || (numel(x)>=1)));
    p.parse(varargin{:});

    updateRTVFobjectAtPosition = p.Results.updateRTVFobjectAtPosition;
    updateRTVFobjectWithCenterConesNum = p.Results.updateRTVFobjectWithCenterConesNum;

    % Generate mosaic filename and directory
    [mosaicFileName, mosaicDirectory] = midgetRGCMosaicInspector.generateMosaicFileName(mosaicCenterParams);

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
        [~,updatedRTVFobjectIndex] = min(d);
        eccentricitySamplingGrid = eccentricitySamplingGrid(updatedRTVFobjectIndex,:);
        fprintf(2, 'Will update the RTVFobjList at index %d (pos (degs) = %2.2f,%2.2f)', ...
            updatedRTVFobjectIndex, eccentricitySamplingGrid(1), eccentricitySamplingGrid(2));
    end


    % Visualize the mosaic
    theMidgetRGCmosaic.visualize( ...
        'eccentricitySamplingGrid', eccentricitySamplingGrid, ...
        'inputPoolingVisualization', 'centerOnly');

    % Assemble R2CFT filename
    R2VFTobjFileName = midgetRGCMosaicInspector.R2VFTobjFileNameForMosaicFileName(mosaicFileName, mosaicSurroundParams.H1cellIndex);

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
        % Generate list of RTVT objects
        [theRTFVTobjList, theOpticsPositionGrid, ...
         theConesNumPooledByTheRFcenterGrid, ...
         theVisualSTFSurroundToCenterRcRatioGrid, ...
         theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid] = midgetRGCMosaic.R2VFTobjects(...
                    theMidgetRGCmosaic, eccentricitySamplingGrid, ...
                    mosaicSurroundParams, opticsParams, fitParams);
    else
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

        for iUpdatedIndex = 1:numel(updateRTVFobjectWithCenterConesNum)

            % Compute sourceRTVFobjectIndex
            sourceRTVFobjectIndex = find(theUpdatedConesNumPooledByTheRFcenterGrid == updateRTVFobjectWithCenterConesNum(iUpdatedIndex));

            % Compute destinationRTVFobjectIndex
            % Retrieve the indices of the fitted RTVF objects that have the
            % same # of center cones
            centerConeMatchObjIndices = find(theConesNumPooledByTheRFcenterGrid == updateRTVFobjectWithCenterConesNum(iUpdatedIndex));

            % Compute distance based weights for this RGC and the fitted RTVF objects
            distancesToSamplingGridPositions = sqrt(sum((bsxfun(@minus, theOpticsPositionGrid(centerConeMatchObjIndices,:), targetPosition)).^2,2));
            [~, idx] = min(distancesToSamplingGridPositions);
            destinationRTVFobjectIndex = centerConeMatchObjIndices(idx);

            % Query user whether to updaste the list of RVFTobj
            prompt = sprintf('Update the RTVFlist with current RTVFobj with %d cones at position (degs): %2.2f %2.2f? [y/n] : ', ...
                theConesNumPooledByTheRFcenterGrid(destinationRTVFobjectIndex), ...
                theOpticsPositionGrid(destinationRTVFobjectIndex,1), theOpticsPositionGrid(destinationRTVFobjectIndex,2));
            txt = lower(input(prompt,'s'));
            if isempty(txt)
                txt = 'n';
            end
            if (strcmp(txt, 'n'))
                fprintf('Will skip overwriting previous data.\n');
                continue;
            end
            fprintf('Overwriting previous data.\n');

            % Update !
            theRTFVTobjList{destinationRTVFobjectIndex} = theUpdatedRTFVTobjList{sourceRTVFobjectIndex};
        end

        % Save the updated RTVFT list
        save(R2VFTobjFileName, ...
            'theRTFVTobjList', ...
            'theOpticsPositionGrid', 'theConesNumPooledByTheRFcenterGrid', ...
            'theVisualSTFSurroundToCenterRcRatioGrid', ...
            'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid');
    else
        % Save the computed RTVFT list
        save(R2VFTobjFileName, 'theRTFVTobjList', 'theOpticsPositionGrid', ...
                    'theConesNumPooledByTheRFcenterGrid', ...
                    'theVisualSTFSurroundToCenterRcRatioGrid', ...
                    'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid');
        fprintf('Computed R2VFTobjects saved in: %s\n', R2VFTobjFileName);
    end
end