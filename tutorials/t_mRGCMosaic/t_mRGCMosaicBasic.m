%% Introduction to the midget RGC mosaic (mRGCMosaic) object.
%
% Description:
%    Demonstrates
%        - creation o a midget RGC mosaic,
%        - how to compute with it, and
%        - how to visualize different aspects of the mRGCMosaic
%        - how to visualize its response
%


% History:
%    01/27/23  NPC  ISETBIO Team, Copyright 2023 Wrote it.

function t_mRGCMosaicBasic

%% Load the source midgetRGCMosaic from the database of prebaked mRGCMosaics
% Choose the prebaked mRGCMosaic that was generated at an eccentricity of (0,0) 
% extending over a 3x3 deg, with cone pooling weights tuned for the Polans subject
% with rank 6 and a pupil size of 3.0 mm

mosaicCenterParams = struct(...
    'positionDegs',[0 0], ...
    'sizeDegs',  [3 3], ...        
    'whichEye', 'right eye');

opticsParams = struct(...
            'ZernikeDataBase', 'Polans2015', ...
            'subjectRankOrder', 6, ...
            'pupilDiameterMM', 3.0 ...
        );

rfModelParams = struct(...
    'H1cellIndex', 1 ...
    );

% Load it
theSourceMidgetRGCMosaic = loadSourceMidgetRGCMosaic(mosaicCenterParams, rfModelParams, opticsParams);

% The SourceMidgetRGCMosaic contains all the information that was used to
% derive the weights. It is not compute-ready.

% Instantiate a compute-ready mRGCMosaic from the sourceMidgetRGCMosaic
% Here we are using part of the sourceMidgetRGCMosaic, centered at (x,y) = (1,0.5), 
% with width = 0.4 degs and height = 0.2 degs

theMRGCMosaic = mRGCMosaic(theSourceMidgetRGCMosaic, ...
        'eccentricityDegs', [1 0.5], ...
        'sizeDegs', [0.5 0.5], ...
        'name', 'my small off-center mRGC mosaic', ...
        'beVerbose', true, ...
        'visualizeSpatialRelationshipToSourceMosaic', true);


theMRGCMosaic2 = mRGCMosaic(theSourceMidgetRGCMosaic, ...
        'eccentricityDegs', [1 1], ...
        'sizeDegs', [2 2], ...
        'name', 'my large mRGC mosaic', ...
        'beVerbose', true, ...
        'visualizeSpatialRelationshipToSourceMosaic', true);



oneCenterConeRGCindices = find(theMRGCMosaic2.centerSubregionConesNums == 1);
twoCenterConeRGCindices = find(theMRGCMosaic2.centerSubregionConesNums == 2);
visualizedRGCindices = [oneCenterConeRGCindices(1:2) twoCenterConeRGCindices(1:2)]

% Visualize the centers of all RGCs, identifying 2 1-cone center RGCs and 2, 2-cone center RGCs
theMRGCMosaic2.visualize(...,
    'identifyInputCones', false, ...
    'labelRGCsWithIndices', visualizedRGCindices, ...
    'backgroundColor', [0.7 0.7 0.7]);


theMRGCMosaic2.multifocalRTVFgrids
theMRGCMosaic2.multifocalRTVFopticsParams

% Visualize the RFs of the identified  RGCs
theMRGCMosaic2.visualizeRFs(visualizedRGCindices);

end


% ======== HELPER FUNCTIONS ========

function theSourceMidgetRGCMosaic = loadSourceMidgetRGCMosaic(mosaicCenterParams, rfModelParams, opticsParams)

    dropboxDir = midgetRGCMosaicInspector.localDropboxPath;
    frozenMidgetRGCMosaicsDir = 'productionMidgetRGCMosaics/frozenMosaics';
    directoryPath = fullfile(dropboxDir,frozenMidgetRGCMosaicsDir);

    
    sourceMidgetRGCMosaicFileName = ...
        sprintf('MRGCmosaic_Ecc_%2.1f_%2.1f_sizeDegs_%2.1f_%2.1f_%s_Rank_%d_Pupil_%2.1f_H1cellIndex_%d_Frozen.mat', ...
        mosaicCenterParams.positionDegs(1), mosaicCenterParams.positionDegs(2), ...
        mosaicCenterParams.sizeDegs(1), mosaicCenterParams.sizeDegs(2), ...
        opticsParams.ZernikeDataBase, ...
        opticsParams.subjectRankOrder, ...
        opticsParams.pupilDiameterMM, ...
        rfModelParams.H1cellIndex);

    theFileName = fullfile(directoryPath, sourceMidgetRGCMosaicFileName);
    fprintf('Will try to load %s ... \n', sourceMidgetRGCMosaicFileName)
    fprintf('from %s ... \n', directoryPath);

    % Check that the mosaic directory exists
    assert(isfolder(directoryPath), sprintf('Mosaic directory (''%s'') not found.', directoryPath));

    % Check that the mosaic file exists
    assert(isfile(theFileName), sprintf('Mosaic file (''%s'') not found.', theFileName));

    % Mosaic file found, so load the data
    load(theFileName, 'theMidgetRGCmosaic');
    theSourceMidgetRGCMosaic = theMidgetRGCmosaic;
    clear 'theMidgetRGCmosaic';
    fprintf('Loaded source mosaic.\n');
end
