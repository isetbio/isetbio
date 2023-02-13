function generateR2VFTobjects(mosaicCenterParams, rfModelParams, opticsParams, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('RTVobjIndicesToBeComputed', 'all',  @(x)(((ischar(x))&&(strcmp(x, 'all'))) || isnumeric(x)));
    p.parse(varargin{:});

    RTVobjIndicesToBeComputed = p.Results.RTVobjIndicesToBeComputed;

    if (strcmp(RTVobjIndicesToBeComputed, 'all'))
        midgetRGCMosaicInspector.say('Generating select RTVF objects');
    else
        midgetRGCMosaicInspector.say('Generating all RTVF objects');
    end

    % Generate mosaic filename and directory
    [mosaicFileName, mosaicDirectory] = midgetRGCMosaicInspector.mosaicFileName(mosaicCenterParams);

    fprintf('Loading the midget RGC mosaic from %s.\nPlease wait ...\n', mosaicFileName);
    % Load the center-connected mosaic
    load(mosaicFileName, 'theMidgetRGCmosaic');

    % Generate the eccentricitySamplingGrid
    eccentricitySamplingGrid = midgetRGCMosaic.eccentricitySamplingGridCoords(...
        mosaicCenterParams.positionDegs, mosaicCenterParams.sizeDegs, ...
        theMidgetRGCmosaic.rgcRFpositionsDegs, ...
        rfModelParams.eccentricitySamplingGridHalfSamplesNum, ...
        'hexagonal', true);

    fprintf('\nVisualizing the mosaic. Please wait ...')
    % Visualize the mosaic
    theMidgetRGCmosaic.visualize( ...
        'eccentricitySamplingGrid', eccentricitySamplingGrid, ...
        'maxVisualizedRFs', 0);
    fprintf('\nDone ! \n');

     % Assemble R2CVFT filename
    R2VFTobjFileName = midgetRGCMosaicInspector.R2VFTobjFileName(mosaicFileName, opticsParams, rfModelParams.H1cellIndex);

    % Ask the user if he wants to use a dictionary with previously
    % fitted params to use as initial values
    usePreviousFittedParamsValues = input('Use initial values from a previous fit ? [y = YES] ', 's');
    initialGridRetinalConePoolingParamsStruct = [];
    if strcmpi(usePreviousFittedParamsValues, 'y')
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
    
    fitParams = struct();
    fitParams.exportsDirectory = mosaicDirectory;
    fitParams.initialGridRetinalConePoolingParamsStruct = initialGridRetinalConePoolingParamsStruct;

    tStart = cputime;

    midgetRGCMosaic.R2VFTobjects(...
                    RTVobjIndicesToBeComputed, ...
                    theMidgetRGCmosaic, eccentricitySamplingGrid, ...
                    rfModelParams, opticsParams, fitParams);

    timeLapsedMinutes = (cputime - tStart)/60;
    fprintf('\n\n midgetRGCMosaic.R2VFTobjects were generated in %d positions and fitting took %f minutes\n', ...
            size(eccentricitySamplingGrid,1), timeLapsedMinutes);

    % Save the computed RTVFT list
    save(R2VFTobjFileName, 'theRTFVTobjList', 'theOpticsPositionGrid', ...
                    'theConesNumPooledByTheRFcenterGrid', ...
                    'theVisualSTFSurroundToCenterRcRatioGrid', ...
                    'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid', ...
                    '-v7.3');
    fprintf('Computed R2VFTobjects saved in: %s\n', R2VFTobjFileName);

end



