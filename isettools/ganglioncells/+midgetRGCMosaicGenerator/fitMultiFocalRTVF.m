function fitMultiFocalRTVF(mosaicCenterParams, rfModelParams, opticsParams, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('RTVobjIndicesToBeComputed', 'all',  @(x)(((ischar(x))&&(strcmp(x, 'all'))) || isnumeric(x)));
    p.addParameter('computeLconeCenterComputeStruct', true, @islogical);
    p.addParameter('computeMconeCenterComputeStruct', true, @islogical);
    p.addParameter('multiStartsNumRetinalPooling', 1, @isscalar);
    p.parse(varargin{:});

    RTVobjIndicesToBeComputed = p.Results.RTVobjIndicesToBeComputed;
    computeLconeCenterComputeStruct = p.Results.computeLconeCenterComputeStruct;
    computeMconeCenterComputeStruct = p.Results.computeMconeCenterComputeStruct;
    multiStartsNumRetinalPooling = p.Results.multiStartsNumRetinalPooling;

    if (ischar(RTVobjIndicesToBeComputed))&&(strcmp(RTVobjIndicesToBeComputed, 'all'))
        midgetRGCMosaicInspector.say('Generating all RTVF objects');
        fullWrite = true;
    else
        midgetRGCMosaicInspector.say('Generating select RTVF objects');
        fullWrite = false;
    end

    % Generate mosaic filename and directory
    [mosaicFileName, mosaicDirectory] = midgetRGCMosaicInspector.mosaicFileName(mosaicCenterParams);

    fprintf('Loading the midget RGC mosaic from %s.\nPlease wait ...\n', mosaicFileName);
    % Load the center-connected mosaic
    load(mosaicFileName, 'theMidgetRGCmosaic');

    % Instantiate the multifocalRTVFobj
    samplingScheme = 'hexagonal';

    % Instantiate an @RTVFmultifocal object with the desired mosaic and
    % optics and rfModel params. During instantiation the @RTVFmultifocal
    % generates a grid of spatial positions/center cones num on which a
    % different @RTVF object will be fit, separately for L- and for
    % M-center RFs
    theMultifocalRTVFOBJ = RTVFmultifocal(theMidgetRGCmosaic, ...
        mosaicCenterParams, opticsParams, rfModelParams, ...
        samplingScheme, ...
        'multiStartsNumRetinalPooling', multiStartsNumRetinalPooling);

    fprintf('\nVisualizing the mosaic. Please wait ...')
    % Visualize the mosaic
    theMidgetRGCmosaic.visualize( ...
        'eccentricitySamplingGrid', theMultifocalRTVFOBJ.nominalSpatialSamplingGrid, ...
        'maxVisualizedRFs', 0);
    fprintf('\nDone ! \n');
    

    % Ask the user if he wants to use a dictionary with previously
    % fitted params to use as initial values
    usePreviousFittedParamsValues = midgetRGCMosaicInspector.queryUserForYesNoResponse('\nUse initial values from a previous fit ?');
    initialGridRetinalConePoolingParamsStruct = [];
    if (usePreviousFittedParamsValues)
        dropboxDir = midgetRGCMosaicInspector.localDropboxPath();

        midgetRGCMosaicInspector.say('Select file with previously derived retinal cone pooling params');

        [file,path] = uigetfile(fullfile(dropboxDir, '*.mat'), ...
                            'Select a file');
    
        if (file ~= 0)
            fName = fullfile(path,file);
            load(fName, 'retinalConePoolingParamsDictionary', 'theConesNumPooledByTheRFcenterGrid', 'theOpticsPositionGrid');
            initialGridRetinalConePoolingParamsStruct.dictionary = retinalConePoolingParamsDictionary;
            initialGridRetinalConePoolingParamsStruct.eccentricitySamplingGrid = theOpticsPositionGrid;
            initialGridRetinalConePoolingParamsStruct.centerConesNumGrid = theConesNumPooledByTheRFcenterGrid;
        end
    end
    
    % Go ! Fit a separate @RTVF object to each node of this grid.
    tStart = cputime;
    theMultifocalRTVFOBJ.compute( ...
        initialGridRetinalConePoolingParamsStruct, ...
        RTVobjIndicesToBeComputed, ...
        computeLconeCenterComputeStruct, ...
        computeMconeCenterComputeStruct, ...
        mosaicDirectory);


    timeLapsedMinutes = (cputime - tStart)/60;
    fprintf('\n\n midgetRGCMosaic.R2VFTobjecs fitting took %f minutes\n', ...
           timeLapsedMinutes);

    
    % Assemble R2CVFT filename
    R2VFTobjFileName = midgetRGCMosaicInspector.R2VFTobjFileName(mosaicFileName, opticsParams, rfModelParams.H1cellIndex);

    if (fullWrite)
        % Save everything
        % Save the computed RTVFT list
        theRTFVTobjList = theMultifocalRTVFOBJ.RTVFTobjList;
        theOpticsPositionGrid = theMultifocalRTVFOBJ.opticalPositionGrid;
        theConesNumPooledByTheRFcenterGrid = theMultifocalRTVFOBJ.conesNumPooledByTheRFcenterGrid;
        theVisualSTFSurroundToCenterRcRatioGrid = theMultifocalRTVFOBJ.visualSTFSurroundToCenterRcRatioGrid;
        theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid = theMultifocalRTVFOBJ.visualSTFSurroundToCenterIntegratedSensitivityRatioGrid;

        save(R2VFTobjFileName, 'theRTFVTobjList', ...
                        'theOpticsPositionGrid', ...
                        'theConesNumPooledByTheRFcenterGrid', ...
                        'theVisualSTFSurroundToCenterRcRatioGrid', ...
                        'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid', ...
                        '-v7.3');
        fprintf('All computed R2VFTobjects were saved in: %s\n', R2VFTobjFileName);
    else
        % Ask the user whether to overwrite the R2VFTobjFileName
        R2VFTobjFileName
        updateFile = midgetRGCMosaicInspector.queryUserForYesNoResponse('\nUpdate the all RTVF objects file with the re-fitted RTVFs?');
        if (updateFile)
            % Load previous file with an existing theRTFVTobjList
            load(R2VFTobjFileName, 'theRTFVTobjList', ...
                        'theOpticsPositionGrid', ...
                        'theConesNumPooledByTheRFcenterGrid', ...
                        'theVisualSTFSurroundToCenterRcRatioGrid', ...
                        'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid');

            for iObj = 1:numel(RTVobjIndicesToBeComputed)
                theOverWrittenOBJindex = RTVobjIndicesToBeComputed(iObj);
                fprintf('Overwrting RTVF obj #%d\n', theOverWrittenOBJindex);
                theRTFVTobjList{theOverWrittenOBJindex} = theMultifocalRTVFOBJ.RTVFTobjList{theOverWrittenOBJindex};
            end

            save(R2VFTobjFileName, 'theRTFVTobjList', ...
                        'theOpticsPositionGrid', ...
                        'theConesNumPooledByTheRFcenterGrid', ...
                        'theVisualSTFSurroundToCenterRcRatioGrid', ...
                        'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid', ...
                        '-v7.3');
            fprintf('Updated the all RTVF objects file (%s)\n', R2VFTobjFileName);

        else
            fprintf('Did not update the all RTVF objects file ');
        end

    end

end