
% Initialize session
close all; clear all;

% Function handle to derive from linear distance in MMs to angular distance
% in degrees on the macaque retina (222 microns/deg)
micronsPerDegreeInMacaqueRetina = 222;
macaqueLinearDistanceMMtoAngularDistanceDegs = @(linearMMs)(linearMMs*1e3/ micronsPerDegreeInMacaqueRetina);


% Load center-connected MRGC mosaic at 28 degs
mosaicEccDegs = [-28 0];
mosaicSizeDegs = 20*[1 1];

extraSizeDegs = 1;
whichEye = 'right eye';

centerConnectedParamsStruct.whichEye = whichEye;
centerConnectedParamsStruct.eccentricityDegs = mosaicEccDegs;
centerConnectedParamsStruct.sizeDegs = mosaicSizeDegs+extraSizeDegs*[1 1];
centerConnectedParamsStruct.rfCenterConnectivityParams.chromaticSpatialVarianceTradeoff = 1.0;

[theCenterConnectedMRGCMosaicFullFileName, theIntermediateConnectivityStageMetaDataFile, theCenterConnectedMRGCMosaicFileName] = ...
        RGCMosaicConstructor.exportedMosaicFileName(centerConnectedParamsStruct, 'center connected');


load(theCenterConnectedMRGCMosaicFullFileName, 'theMRGCMosaic');


% Target in linear ecc (3.5 mm, temporal retina)
% EJ (2023) paper, Fig 1A
targetLinearEccMM(1) = -3.5;
targetLinearEccMM(2) = 0.0;

% EJ (2023) paper, Fig 1B (8.5 mm TEE)
targetLinearEccMM(1) = -8.4;
targetLinearEccMM(2) = sqrt(8.5^2 - targetLinearEccMM(1)^2);

% Compute the corresponding angular eccentricity for human retina 
targetAngularEccDegsInHumanRetina = theMRGCMosaic.linearEccMMsToAngularEccDegs(targetLinearEccMM)

% Compute the corresponding angular eccentricity for macaque retina 
targetAngularEccDegsInMacaqueRetina = macaqueLinearDistanceMMtoAngularDistanceDegs(targetLinearEccMM)
pause
% Scale bar for RFs specified in linear size (microns)
receptiveFieldScaleBarInMicrons = 21;

% Compute the corresponding angular size for human retina
receptiveFieldScaleBarInDegsInHumanRetina = theMRGCMosaic.linearEccMMsToAngularEccDegs(receptiveFieldScaleBarInMicrons*1e-3);
receptiveFieldScaleBarInDegsInMacaqueRetina = macaqueLinearDistanceMMtoAngularDistanceDegs(receptiveFieldScaleBarInMicrons*1e-3);

% Scale bar for mosaics specified in linear size (microns)
mosaicScaleBarInMicrons = 500;

% Compute the corresponding angular size for human retina
mosaicScaleBarInDegsInHumanRetina = theMRGCMosaic.linearEccMMsToAngularEccDegs(mosaicScaleBarInMicrons*1e-3);
mosaicScaleBarInDegsInMacaqueRetina = macaqueLinearDistanceMMtoAngularDistanceDegs(mosaicScaleBarInMicrons*1e-3);


%scenario = 'mosaic @ angular ecc == to macaque angular ecc of the target position (retinal mm)';
scenario = 'mosaic @ angular ecc corresponding to the target position (retinal mm)';

switch (scenario)
    case 'mosaic @ angular ecc == to macaque angular ecc of the target position (retinal mm)'
        % Depict the human mosaic at an angular eccentricity that is equal to the 
        % **macaque retina** angular eccentricity that corresponds to the targetLinearEccMM
        xo = targetAngularEccDegsInMacaqueRetina(1)
        yo = targetAngularEccDegsInMacaqueRetina(2)
        scaleBarDegsSingleRF = receptiveFieldScaleBarInDegsInHumanRetina;
        scaleBarDegsMosaic = mosaicScaleBarInDegsInHumanRetina;

    case 'mosaic @ angular ecc corresponding to the target position (retinal mm)'
        % Depict the human mosaic at an angular eccentricity that corresponds 
        % to the targetLinearEccMM
        xo = targetAngularEccDegsInHumanRetina(1)
        yo = targetAngularEccDegsInHumanRetina(2)
        scaleBarDegsSingleRF = receptiveFieldScaleBarInDegsInHumanRetina;
        scaleBarDegsMosaic = mosaicScaleBarInDegsInHumanRetina;

    otherwise
        error('Unknown scenario: ''%s''.', scenario);
end


domainVisualizationLimits = [xo-2.0 xo+2.0 yo-2.0 yo+2.0];
domainVisualizationTicks = struct('x', -50:0.5:50, 'y', -50:0.5:50);

d2 = sum((bsxfun(@minus, theMRGCMosaic.rgcRFpositionsDegs, [xo yo])).^2,2);
[~,theRGCindex] = min(d2(:));

hFig = figure(1000); clf;
% Prepare figure and axes
ff = PublicationReadyPlotLib.figureComponents('1x1 giant square mosaic');
theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);

theMRGCMosaic.visualize(...
    'figureHandle', hFig, ...
    'axesHandle', theAxes{1,1}, ...
    'scaleBarDegs', scaleBarDegsMosaic, ...
    'labelRGCsWithIndices', theRGCindex, ...
    'labeledRGCsColor', [1 0 0], ...
    'domainVisualizationLimits', domainVisualizationLimits, ...
    'domainVisualizationTicks', domainVisualizationTicks, ...
    'plottedRFoutlineLineWidth', 1.0, ...
    'plotTitle', scenario);

% Finalize figure using the Publication-Ready format
PublicationReadyPlotLib.applyFormat(theAxes{1,1},ff);
    
% Export figure
theRawFiguresDir = RGCMosaicConstructor.rawFigurePDFsDir();
thePDFfileName = fullfile(theRawFiguresDir, strrep(theCenterConnectedMRGCMosaicFileName, '.mat', sprintf('roi_%2.1fdegs_%s.pdf',xo, scenario) ));
NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);


hFig = figure(1001); clf;
% Prepare figure and axes
ff = PublicationReadyPlotLib.figureComponents('1x1 giant square mosaic');
theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);


% Prepare figure and axes
ff = PublicationReadyPlotLib.figureComponents('1x1 standard figure');
theAxes = PublicationReadyPlotLib.generatePanelAxes(hFig,ff);


theMRGCMosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', theAxes{1,1}, ...
            'scaleBarDegs', scaleBarDegsSingleRF, ...
            'singleRGCconePoolingRFmaps', true, ...
            'visualizedRGCindices', theRGCindex, ...
            'tickSeparationArcMinForRFconePoolingMap', 5, ...
            'plotTitle', scenario)

% Finalize figure using the Publication-Ready format
PublicationReadyPlotLib.applyFormat(theAxes{1,1},ff);
    
% Export figure
theRawFiguresDir = RGCMosaicConstructor.rawFigurePDFsDir();
thePDFfileName = fullfile(theRawFiguresDir, strrep(theCenterConnectedMRGCMosaicFileName, '.mat', sprintf('roi_%2.1fdegs_RGC%d_RF_%s.pdf',xo, theRGCindex, scenario)));
NicePlot.exportFigToPDF(thePDFfileName,hFig,  300);
