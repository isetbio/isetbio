%
% pStruct = RGCMosaicConstructor.helper.utils.initializeRGCMosaicGenerationParameters(...
%	coneMosaicSpecies, opticsSubjectName);
%

function pStruct = initializeRGCMosaicGenerationParameters(coneMosaicSpecies, ...
    opticsSubjectName, rgcMosaicName, targetVisualSTFdescriptor)

% ----- CONE MOSAIC PARAMS ----

% Which mRGC/cone mosaic lattice to use
pStruct.sourceLatticeSizeDegs = 64; 
pStruct.whichEye = 'right eye';

switch (coneMosaicSpecies)
	case 'human'
		% Human retina
		pStruct.customLMSconeDensities = [];

	case 'macaque'
		% Macaque retina, 1:1 L:M cone ratio, 10% S-cones
		pStruct.customLMSconeDensities = [0.45 0.45 0.1];

	otherwise
		error('Unknown cone mosaic species: ''%s''.', coneMosaicSpecies);
end % switch


% ----- OPTICS PARAMS ----

pStruct.whichZernikeDataBase = 'Polans2015';

if (isscalar(opticsSubjectName))
	pStruct.whichSubjectID = opticsSubjectName;
elseif (ischar(opticsSubjectName))
	switch (opticsSubjectName)
		case {'PLOSpaperDefaultSubject', 'PLOSpaperStrehlRatio_0.59'}
			% This is ranked #3 in terms of resolution at the fovea (45 c/deg), Strehl Ratio around 0.59
			pStruct.whichSubjectID = 2;

		case 'PLOSpaperStrehlRatio_0.87'
			% This is ranked #1 in terms of resolution at the fovea (61 c/deg), Strehl Ratio around 0.87
			pStruct.whichSubjectID = 6;

		case 'PLOSpaperStrehlRatio_0.72'
			% This is ranked #6 in terms of resolution at the fovea (31 c/deg), Strehl Ratio around 0.72
			pStruct.whichSubjectID = 4;
	
		case 'PLOSpaperStrehlRatio_0.60'
			% This is ranked #4 in terms of resolution at the fovea (36 c/deg),  Strehl Ratio around 0.60
			pStruct.whichSubjectID = 10;

		case 'PLOSpaperStrehlRatio_0.27'
			% This is ranked #5 in terms of resolution at the fovea (31 c/deg), Strehl Ratio around 0.27
			pStruct.whichSubjectID = 8;

		case 'PLOSpaperStrehlRatio_0.23'
			% This is ranked #9 in terms of resolution at the fovea (18 c/deg), Strehl Ratio around 0.23
			pStruct.whichSubjectID = 1;

        case {'PLOSpaperSecondSubject', 'PLOSpaperStrehlRatio_0.21'}
			% This is ranked #7 in terms of resolution at the fovea (23 c/deg),  Strehl Ratio around 0.21
			pStruct.whichSubjectID = 3;

		case 'PLOSpaperStrehlRatio_0.19'
			% This is ranked #8 in terms of resolution at the fovea (19 c/deg), Strehl Ratio around 0.19
			pStruct.whichSubjectID = 5;

		case 'PLOSpaperStrehlRatio_0.09'
			% This is ranked #10 in terms of resolution at the fovea (9 c/deg), Strehl Ratio around 0.09
			pStruct.whichSubjectID = 7;

		case 'VSS2024TalkFirstSubject'
			% This is ranked #7 in terms of resolution at the fovea (23 c/deg)
			pStruct.whichSubjectID = 3;

		case 'VSS2024TalkSecondSubject'
			% This is ranked #5 in terms of resolution at the fovea (31 c/deg)
			pStruct.whichSubjectID = 8;
		otherwise
			error('Unknown optics subject name: ''%s''.', opticsSubjectName);
	end % switch (opticsSubjectName)
else
	opticsSubjectName
	error('Passed opticsSubjectName is neither a scalar nor a string');
end


% ----- RGC MOSAIC PARAMS -----
% RF center spatialChromaticUniformityTradeoff [0: minimize chromatic variance 1: minimize spatial variance]
pStruct.rgcMosaic.spatialChromaticUniformityTradeoff = 1.0;

% Use mosaic with RFcenter overlap
pStruct.rgcMosaic.employRFCenterOverlappingMosaic = true;

% RF surround optimization params
switch (rgcMosaicName)
	case {'VSS2024TalkTemporal3DegsMosaic'}
		rgcMosaicSurroundOptimization.mosaicEccDegs = [-3 0];
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 2*[1 1];
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'high res';
		rgcMosaicSurroundOptimization.minGridSize = 0.5;
		rgcMosaicSurroundOptimization.maxGridSize = 1.5;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.08^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'LowerLeftQH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = false;

    case {'VSS2024TalkTemporal7DegsMosaic'}
        rgcMosaicSurroundOptimization.mosaicEccDegs = [-7 0];
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 4*[1 1];
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'medium res'; 
		rgcMosaicSurroundOptimization.minGridSize = 1.0; 
		rgcMosaicSurroundOptimization.maxGridSize = 3.0;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.16^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'LowH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = false;

	case 'PLOSpaperFovealMosaic'
		rgcMosaicSurroundOptimization.mosaicEccDegs = [0 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = [2 2];    % [1 1] (28182 cells)
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'high res'; 
		rgcMosaicSurroundOptimization.minGridSize = 0.25; 
		rgcMosaicSurroundOptimization.maxGridSize = 1.0;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.05^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'LowH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = false;

	case 'PLOSpaperTemporal2DegsMosaic' 
		rgcMosaicSurroundOptimization.mosaicEccDegs = [-2 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 2*[1 1]; % [-1 -3] (10366 cells)
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'high res'; 
		rgcMosaicSurroundOptimization.minGridSize = 0.5; 
		rgcMosaicSurroundOptimization.maxGridSize = 1.5;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.05^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'LowH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = false;

	case 'PLOSpaperTemporal4DegsMosaic'
		rgcMosaicSurroundOptimization.mosaicEccDegs = [-4 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 3*[1 1]; % [-2.5 -5.5] (9303 cells)
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'high res'; 
		rgcMosaicSurroundOptimization.minGridSize = 0.5; 
		rgcMosaicSurroundOptimization.maxGridSize = 1.5;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.05^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'LowH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = false;


	case 'PLOSpaperTemporal7DegsMosaic'
		rgcMosaicSurroundOptimization.mosaicEccDegs = [-7 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 4*[1 1]; %-5 -9 (6979 cells)
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'medium res'; 
		rgcMosaicSurroundOptimization.minGridSize = 1.0; 
		rgcMosaicSurroundOptimization.maxGridSize = 3.0;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.1^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'LowH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = false;

	case 'PLOSpaperTemporal10DegsMosaic'
		rgcMosaicSurroundOptimization.mosaicEccDegs = [-10 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 4*[1 1]; % -8 -12 (3772 cells) (linear part of super Gaussian exponent)
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'medium res'; 
		rgcMosaicSurroundOptimization.minGridSize = 1.0; 
		rgcMosaicSurroundOptimization.maxGridSize = 3.0;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.17^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'LowH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = false;

	case 'PLOSpaperTemporal14DegsMosaic'
		rgcMosaicSurroundOptimization.mosaicEccDegs = [-14 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 5*[1 1]; %-11.5 -16.5 (2695 cells)
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'medium res'; 
		rgcMosaicSurroundOptimization.minGridSize = 1.0; 
		rgcMosaicSurroundOptimization.maxGridSize = 3.0;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.18^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'LowH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = false;


	case 'PLOSpaperTemporal19DegsMosaic'
		rgcMosaicSurroundOptimization.mosaicEccDegs = [-19 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 6*[1 1];
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'medium res'; 
		rgcMosaicSurroundOptimization.minGridSize = 1.5; 
		rgcMosaicSurroundOptimization.maxGridSize = 3.0;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.18^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'MidH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = true;

	case 'PLOSpaperTemporal25DegsMosaic'
		rgcMosaicSurroundOptimization.mosaicEccDegs = [-25 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 7*[1 1]; %-21.5 -28.5 (1358 cells)
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'low res'; 
		rgcMosaicSurroundOptimization.minGridSize = 2; 
		rgcMosaicSurroundOptimization.maxGridSize = 4.0;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.15^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'UpperH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = true;

    case 'PLOSpaperTemporal32DegsMosaic'
		rgcMosaicSurroundOptimization.mosaicEccDegs = [-32 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 9*[1 1]; %-27.5 -36.5 (1358 cells)
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'low res'; 
		rgcMosaicSurroundOptimization.minGridSize = 2; 
		rgcMosaicSurroundOptimization.maxGridSize = 4.0;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.15^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'UpperH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = true;


    case 'PLOSpaperNasal2DegsTinyMosaic'
        rgcMosaicSurroundOptimization.mosaicEccDegs = [2 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 0.5*[1 1]; 
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'medium res'; 
        rgcMosaicSurroundOptimization.minGridSize = 0.5; 
		rgcMosaicSurroundOptimization.maxGridSize = 1.5;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.05^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'LowH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = false;


    case 'PLOSpaperNasal7DegsMosaic'
		rgcMosaicSurroundOptimization.mosaicEccDegs = [7 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 11*[1 1]; 
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'medium res'; 
		rgcMosaicSurroundOptimization.minGridSize = 1.0; 
		rgcMosaicSurroundOptimization.maxGridSize = 3.0;
		rgcMosaicSurroundOptimization.addEightExtremePositions = true;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.1^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'LowH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = false;


    case 'PLOSpaperNasal10DegsMosaic'
		rgcMosaicSurroundOptimization.mosaicEccDegs = [10 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 11*[1 1];
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'medium res'; 
		rgcMosaicSurroundOptimization.minGridSize = 1.0; 
		rgcMosaicSurroundOptimization.maxGridSize = 3.0;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.17^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'LowH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = false;

    case 'PLOSpaperNasal14DegsMosaic'
        % This includes the OD, so we have to generate a large enough region
		rgcMosaicSurroundOptimization.mosaicEccDegs = [14 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 13*[1 1]; 
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'medium res'; 
		rgcMosaicSurroundOptimization.minGridSize = 1.0; 
		rgcMosaicSurroundOptimization.maxGridSize = 3.0;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.18^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'LowH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = false;


    case 'PLOSpaperNasal19DegsMosaic'
        % This includes the OD, so we have to generate a large enough region
		rgcMosaicSurroundOptimization.mosaicEccDegs = [19 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 17*[1 1];
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'medium res'; 
		rgcMosaicSurroundOptimization.minGridSize = 1.5; 
		rgcMosaicSurroundOptimization.maxGridSize = 3.0;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.18^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'MidH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = true;


    case 'PLOSpaperNasal25DegsMosaic'
		rgcMosaicSurroundOptimization.mosaicEccDegs = [25 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 7*[1 1]; 
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'low res'; 
		rgcMosaicSurroundOptimization.minGridSize = 2; 
		rgcMosaicSurroundOptimization.maxGridSize = 4.0;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.15^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'UpperH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = true;



	otherwise
		error('No params specified for rgc mosaic  name: ''%s''.', rgcMosaicName);
end % switch (rgcMosaicName)


% Add the specific targetVisualSTFdescriptor
if (isempty(targetVisualSTFdescriptor))
	rgcMosaicSurroundOptimization.targetVisualSTFdescriptor = 'default';
else
	rgcMosaicSurroundOptimization.targetVisualSTFdescriptor = targetVisualSTFdescriptor;
end


% Add all the specified surround optimization params
pStruct.rgcMosaicSurroundOptimization = rgcMosaicSurroundOptimization;

end