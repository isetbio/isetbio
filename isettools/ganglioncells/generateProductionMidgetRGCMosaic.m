function generateProductionMidgetRGCMosaic()

    mosaicCenterParams = struct(...
        'positionDegs',[0 0], ...
        'sizeDegs',  [3 3], ...        
        'whichEye', 'right eye');


    % Generate mosaic filename and directory
    [mosaicFileName, mosaicDirectory] = generateMosaicFileName(mosaicCenterParams);

    % Actions to perform
    actionToPerform = 'generateCenterConnectedMosaic';
    actionToPerform = 'generateR2VFTobjects';

    switch (actionToPerform)
        case 'generateCenterConnectedMosaic'
            generateCenterConnectedMosaic(mosaicCenterParams, mosaicFileName);

        case 'generateR2VFTobjects'

            H1cellIndex = 1;
            mosaicSurroundParams = struct(...
                'eccentricitySamplingGridHalfSamplesNum', 1, ...                         % generate R2VFTobjects at 2*gridHalfSamplesNum + 1 spatial positions
                'centerConnectableConeTypes', [cMosaic.LCONE_ID cMosaic.MCONE_ID], ...   % cone types that can connect to the RF center
                'surroundConnectableConeTypes', [cMosaic.LCONE_ID cMosaic.MCONE_ID], ... % cone types that can connect to the RF surround
                'coneWeightsCompensateForVariationsInConeEfficiency', true, ...          % cone weight compensationfor eccentricity-dependent variations in cone efficiency
                'visualRFmodel', 'gaussian center, gaussian surround', ...
                'retinalConePoolingModel', sprintf('arbitrary center cone weights, double exponential surround from H1 cell with index %d', H1cellIndex),...
                'H1cellIndex', H1cellIndex, ...
                'targetSTFmatchMode', 'STFDoGparams' ...
            );
        
            opticsParams = struct(...
                'ZernikeDataBase', 'Polans2015', ...
                'subjectRankOrder', 6, ...
                'pupilDiameterMM', 3.0 ...
            );

            generateR2VFTobjects(mosaicCenterParams, mosaicSurroundParams, opticsParams, mosaicFileName, mosaicDirectory);

        otherwise
            error('Unknown action: ''%s''.', actionToPerform);
    end


end

function generateR2VFTobjects(mosaicCenterParams, mosaicSurroundParams, opticsParams, mosaicFileName, mosaicDirectory)

    % Load the center-connected mosaic
    load(mosaicFileName, 'theMidgetRGCmosaic');

    
    % Generate the eccentricitySamplingGrid
    eccentricitySamplingGrid = midgetRGCMosaic.eccentricitySamplingGridCoords(...
        mosaicCenterParams.positionDegs, mosaicCenterParams.sizeDegs, ...
        mosaicSurroundParams.eccentricitySamplingGridHalfSamplesNum, ...
        'hexagonal', true);


    %centerConesNum = full(sum(theMidgetRGCmosaic.rgcRFcenterConeConnectivityMatrix,1))
    %[max(centerConesNum) min(centerConesNum) numel(find(centerConesNum ==max(centerConesNum))) numel(find(centerConesNum ==min(centerConesNum)))]

    % Visualize the mosaic
    theMidgetRGCmosaic.visualize( ...
        'eccentricitySamplingGrid', eccentricitySamplingGrid, ...
        'inputPoolingVisualization', 'centerOnly');


    % multiStartsNum: select from:
    % - 1 (Single start run, fastest results), 
    % - some number (Multi-start), or 
    % - inf (Global search)
    fitParams.multiStartsNumDoGFit = 128;

    % Where to save the fitted RTVFobjects
    fitParams.exportsDirectory = mosaicDirectory;

    tStart = cputime;

    % Generate list of RTVT objects
    [theRTFVTobjList, theOpticsPositionGrid, ...
     theConesNumPooledByTheRFcenterGrid, ...
     theVisualSTFSurroundToCenterRcRatioGrid, ...
     theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid] = midgetRGCMosaic.R2VFTobjects(...
                theMidgetRGCmosaic, eccentricitySamplingGrid, ...
                mosaicSurroundParams, opticsParams, fitParams);

    timeLapsedMinutes = (cputime - tStart)/60;
    fprintf('\n\n midgetRGCMosaic.R2VFTobjects were generated in %d positions and fitting took %f minutes\n', timeLapsedMinutes);
    

    % Save the computed list of RTVFTobj and the various grids to the mosaic mat file
    R2VFTobjFileName = R2VFTobjFileNameForMosaicFileName(mosaicFileName, mosaicSurroundParams.H1cellIndex);

    save(R2VFTobjFileName, 'theRTFVTobjList', 'theOpticsPositionGrid', ...
                'theConesNumPooledByTheRFcenterGrid', ...
                'theVisualSTFSurroundToCenterRcRatioGrid', ...
                'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid');
    fprintf('Computed R2VFTobjects saved in: %s\n', R2VFTobjFileName);

end

function generateCenterConnectedMosaic(mosaicParams, mosaicFileName)
    % Generate mRGC mosaic
    theMidgetRGCmosaic = midgetRGCMosaic(...
                'sourceLatticeSizeDegs', 60, ...
                'whichEye', mosaicParams.whichEye, ...
                'eccentricityDegs', mosaicParams.positionDegs, ...
                'sizeDegs', mosaicParams.sizeDegs ...
                );
    % Save the center-connected mosaic
    save(mosaicFileName, 'theMidgetRGCmosaic', '-v7.3');
end

function R2VFTobjFileName = R2VFTobjFileNameForMosaicFileName(mosaicFileName, H1cellIndex)
    R2VFTobjFileName = strrep(mosaicFileName, '.mat', sprintf('_R2VFTobjectsForH1cellIndex_%d.mat', H1cellIndex));
end

function [mosaicFileName, mosaicDirectoryPath] = generateMosaicFileName(mosaicParams)
    dropboxDir = localDropboxPath();
    mosaicDirectoryPath = sprintf('%s/productionMidgetRGCMosaics', dropboxDir);
    mosaicFileName = sprintf('MRGCmosaic_Ecc_%2.1f_%2.1f_sizeDegs_%2.1f_%2.1f.mat', ...
        mosaicParams.positionDegs(1), mosaicParams.positionDegs(2), ...
        mosaicParams.sizeDegs(1), mosaicParams.sizeDegs(2));

    mosaicFileName = fullfile(mosaicDirectoryPath, mosaicFileName);
end


function dropboxDir = localDropboxPath()
    % Get dropboxDir & intermediate data files location
    computerInfo = GetComputerInfo();
    switch (computerInfo.localHostName)
        case 'Ithaka'
            dropboxDir = '/Volumes/SSDdisk/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris';

        case 'Crete'
            dropboxDir = '/Volumes/Dropbox/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris';

        otherwise
            if (contains(computerInfo.networkName, 'leviathan'))
                dropboxDir = '/media/dropbox_disk/Aguirre-Brainard Lab Dropbox/isetbio isetbio';
            else
                error('Could not establish dropbox location')
            end
    end
end
