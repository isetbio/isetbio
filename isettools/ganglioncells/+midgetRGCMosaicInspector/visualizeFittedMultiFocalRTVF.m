function visualizeFittedMultiFocalRTVF(mosaicCenterParams, rfModelParams, opticsParams, varargin)
    % Parse input
    p = inputParser;
    p.addParameter('component', 'all', @(x)(ismember(x, {'fitted locations', 'psfs', 'retinal RFs', 'visual RFs'})));
    p.parse(varargin{:});

    visualizedComponent = p.Results.component;

    % Generate mosaic filename and directory
    [mosaicFileName, mosaicDirectory] = midgetRGCMosaicInspector.mosaicFileName(mosaicCenterParams);

    % Assemble R2CVFT filename
    R2VFTobjFileName = midgetRGCMosaicInspector.R2VFTobjFileName(mosaicFileName, opticsParams, rfModelParams.H1cellIndex);

    % Load the midgetRGCMosaic
    load(mosaicFileName, 'theMidgetRGCmosaic');

    % Load the list of fittted RTVF objects and the sampling grid
    load(R2VFTobjFileName, ...
        'theRTFVTobjList', ...
        'theOpticsPositionGrid', ...
        'theConesNumPooledByTheRFcenterGrid', ...
        'theVisualSTFSurroundToCenterRcRatioGrid', ...
        'theVisualSTFSurroundToCenterIntegratedSensitivityRatioGrid');

    switch (visualizedComponent)
        case 'all'
            figNo = 100;
            RTVFmultifocal.visualizeFittedLocationsCombo(figNo, ...
                theMidgetRGCmosaic, ...
                theRTFVTobjList, ...
                theOpticsPositionGrid, ...
                theConesNumPooledByTheRFcenterGrid, ...
                'bottomPanelInfo', 'retinal M-center RF subregions & M-weighted PSF', ...
                'topPanelInfo', 'retinal L-center RF subregions & L-weighted PSF', ...
                'figPostfix', 'retinalRFs');
        
        case 'fitted locations'
            figNo = 100;
            RTVFmultifocal.visualizeFittedLocations(figNo, ...
                theMidgetRGCmosaic, theOpticsPositionGrid, ...
                theConesNumPooledByTheRFcenterGrid);
        otherwise
            error('Not implemented')
    end

end
