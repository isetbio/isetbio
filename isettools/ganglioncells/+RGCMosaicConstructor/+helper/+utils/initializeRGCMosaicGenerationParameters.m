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

% PSF OBSERVATIONS AT 25 DEG TEMPORAL
% whichSubjectID = 1  % Strehl ratio (25 deg temporal): 0.035 ELONGATED
%whichSubjectID = 2   % Strehl ratio (25 deg temporal): 0.04;ELONGATED
%whichSubjectID = 3   % Strehl ratio (25 deg temporal): 0.04;ELONGATED
%whichSubjectID = 4;  % Strehl ratio (25 deg temporal): 0.045 ELONGATED 
%whichSubjectID = 5;   % Strehl ratio (25 deg temporal): 0.075 ELONGATED
%whichSubjectID = 6;   % Strehl ratio (25 deg temporal): 0.035 ELONGATED
%whichSubjectID = 7;   % Strehl ratio (25 deg temporal): 0.17, NOT ELONGATED - DO THIS
%whichSubjectID = 8;    % Strehl ratio (25 deg temporal): 0.085, NOT ELONGATED - DO THIS
%whichSubjectID = 9;    % Strehl ratio (25 deg temporal): 0.047

% PSF OBSERVATIONS AT 0 DEG
%whichSubjectID = 1;    % Strehl ratio (0 deg temporal): 0.24; NICE FOR BLURRY
%whichSubjectID = 2;    % Strehl ratio (0 deg temporal): 0.42; (THE PLOS 2024) typical subject
%whichSubjectID = 3;    % Strehl ratio (0 deg temporal): 0.2; VERY GOOD FOR BLURRY
%whichSubjectID = 4;    % Strehl ratio (0 deg temporal): 0.75; 
%whichSubjectID = 5;    % Strehl ratio (0 deg temporal): 0.2; ELONGATED
%whichSubjectID = 6;    % Strehl ratio (0 deg temporal): 0.9; VERY SHARP
%whichSubjectID = 7;    % Strehl ratio (0 deg temporal): 0.08; NICE FOR BLURRY
%whichSubjectID = 8;    % Strehl ratio (0 deg temporal): 0.27; NICE FOR BLURRY
%whichSubjectID = 9;    % WEIRD
%whichSubjectID = 10;    % Strehl ratio (0 deg temporal): 0.65;

if (isscalar(opticsSubjectName))
	pStruct.whichSubjectID = opticsSubjectName;
elseif (ischar(opticsSubjectName))
	switch (opticsSubjectName)
		case 'PLOSpaperDefaultSubject'
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

		case 'PLOSpaperStrehlRatio_0.21'
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
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 1*[1 1];
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
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 3*[1 1];
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
		rgcMosaicSurroundOptimization.mosaicSizeDegs = [1 1];    % [1 1] (28182 cells)
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'high res'; 
		rgcMosaicSurroundOptimization.minGridSize = 0.25; 
		rgcMosaicSurroundOptimization.maxGridSize = 1.0;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.05^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'LowH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = false;

    case 'PLOSpaperFovealMosaicLowerOverlap'
		rgcMosaicSurroundOptimization.mosaicEccDegs = [0 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = [1.1 1.1];   
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'high res'; 
		rgcMosaicSurroundOptimization.minGridSize = 0.25; 
		rgcMosaicSurroundOptimization.maxGridSize = 1.0;
		rgcMosaicSurroundOptimization.addEightExtremePositions = true;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.05^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'LowH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = false;


	case 'PLOSpaperTemporal2DegsMosaic' 
		rgcMosaicSurroundOptimization.mosaicEccDegs = [-2 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 1*[1 1]; % [-1 -3] (10366 cells)
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'high res'; 
		rgcMosaicSurroundOptimization.minGridSize = 0.5; 
		rgcMosaicSurroundOptimization.maxGridSize = 1.5;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.05^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'LowH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = false;

    case 'PLOSpaperTemporal2DegsMosaicLowerOverlap' 
		rgcMosaicSurroundOptimization.mosaicEccDegs = [-2.1 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 1*[1 1]; % [-1 -3] (10366 cells)
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
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 2*[1 1]; % [-2.5 -5.5] (9303 cells)
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'high res'; 
		rgcMosaicSurroundOptimization.minGridSize = 0.5; 
		rgcMosaicSurroundOptimization.maxGridSize = 1.5;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.05^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'LowH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = false;

    case 'PLOSpaperTemporal4DegsMosaicLowerOverlap'
		rgcMosaicSurroundOptimization.mosaicEccDegs = [-4.1 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 2*[1 1]; % [-2.5 -5.5] (9303 cells)
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
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 3*[1 1]; %-5 -9 (6979 cells)
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'medium res'; 
		rgcMosaicSurroundOptimization.minGridSize = 1.0; 
		rgcMosaicSurroundOptimization.maxGridSize = 3.0;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.1^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'LowH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = false;

    case 'PLOSpaperTemporal7DegsMosaicLowerOverlap'
		rgcMosaicSurroundOptimization.mosaicEccDegs = [-7.1 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 3*[1 1]; %-5 -9 (6979 cells)
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
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 3*[1 1]; % -8 -12 (3772 cells) (linear part of super Gaussian exponent)
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'medium res'; 
		rgcMosaicSurroundOptimization.minGridSize = 1.0; 
		rgcMosaicSurroundOptimization.maxGridSize = 3.0;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.17^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'LowH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = false;

    case 'PLOSpaperTemporal10DegsMosaicLowerOverlap'
		rgcMosaicSurroundOptimization.mosaicEccDegs = [-10.1 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 3*[1 1]; % -8 -12 (3772 cells) (linear part of super Gaussian exponent)
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
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 4*[1 1]; %-11.5 -16.5 (2695 cells)
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'medium res'; 
		rgcMosaicSurroundOptimization.minGridSize = 1.0; 
		rgcMosaicSurroundOptimization.maxGridSize = 3.0;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.18^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'LowH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = false;

	case 'PLOSpaperTemporal14DegsMosaicLowerOverlap'
		rgcMosaicSurroundOptimization.mosaicEccDegs = [-14.1 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 4*[1 1]; %-11.5 -16.5 (2695 cells)
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
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 5*[1 1];
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'medium res'; 
		rgcMosaicSurroundOptimization.minGridSize = 1.5; 
		rgcMosaicSurroundOptimization.maxGridSize = 3.0;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.18^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'MidH1paramsNarrowVisualSTFparamTolerance'
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = true;

	case 'PLOSpaperTemporal19DegsMosaicLowerOverlap'
		rgcMosaicSurroundOptimization.mosaicEccDegs = [-19.1 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 5*[1 1];
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'medium res'; 
		rgcMosaicSurroundOptimization.minGridSize = 1.5; 
		rgcMosaicSurroundOptimization.maxGridSize = 3.0;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.18^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'MidH1paramsNarrowVisualSTFparamTolerance'
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = true;

	case 'PLOSpaperTemporal25DegsMosaic'
		rgcMosaicSurroundOptimization.mosaicEccDegs = [-25 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 6*[1 1]; %-21.5 -28.5 (1358 cells)
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'low res'; 
		rgcMosaicSurroundOptimization.minGridSize = 2; 
		rgcMosaicSurroundOptimization.maxGridSize = 4.0;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.15^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'UpperH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = true;

    case 'PLOSpaperNasal7DegsMosaic'
		rgcMosaicSurroundOptimization.mosaicEccDegs = [7 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 10*[1 1]; 
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'medium res'; 
		rgcMosaicSurroundOptimization.minGridSize = 1.0; 
		rgcMosaicSurroundOptimization.maxGridSize = 3.0;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.1^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'LowH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = false;


    case 'PLOSpaperNasal10DegsMosaic'
		rgcMosaicSurroundOptimization.mosaicEccDegs = [10 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 10*[1 1];
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
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 12*[1 1]; 
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
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 16*[1 1];
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
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 6*[1 1]; 
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'low res'; 
		rgcMosaicSurroundOptimization.minGridSize = 2; 
		rgcMosaicSurroundOptimization.maxGridSize = 4.0;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.15^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'UpperH1paramsNarrowVisualSTFparamTolerance';
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = true;

	case 'PLOSpaperTemporal25DegsMosaicLowerOverlap'
		rgcMosaicSurroundOptimization.mosaicEccDegs = [-25.1 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 6*[1 1]; %-21.5 -28.5 (1358 cells)
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
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 8*[1 1]; %-27.5 -36.5 (1358 cells)
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'low res'; 
		rgcMosaicSurroundOptimization.minGridSize = 2; 
		rgcMosaicSurroundOptimization.maxGridSize = 4.0;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.15^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'UpperH1paramsNarrowVisualSTFparamTolerance'
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = true;

	case 'PLOSpaperTemporal32DegsMosaicLowerOverlap'
		rgcMosaicSurroundOptimization.mosaicEccDegs = [-32.1 0]; 
		rgcMosaicSurroundOptimization.mosaicSizeDegs = 8*[1 1]; %-27.5 -36.5 (1358 cells)
		rgcMosaicSurroundOptimization.peripheralOptimizationSamplingScheme = 'low res'; 
		rgcMosaicSurroundOptimization.minGridSize = 2; 
		rgcMosaicSurroundOptimization.maxGridSize = 4.0;
		rgcMosaicSurroundOptimization.addEightExtremePositions = false;
		rgcMosaicSurroundOptimization.intSensRatioVariance = 0.15^2; 
		rgcMosaicSurroundOptimization.intSensRatioBias = 1.0;
		rgcMosaicSurroundOptimization.optimizationStrategy = 'UpperH1paramsNarrowVisualSTFparamTolerance'
		rgcMosaicSurroundOptimization.employLconeDominanceOptimizationOnly = true;

	otherwise
		error('Unknown rgc mosaic  name: ''%s''.', rgcMosaicName);
end % switch (rgcMosaicName)


% Extra size (margin)
extraSizeDegs = 1;
rgcMosaicSurroundOptimization.mosaicSizeDegs = rgcMosaicSurroundOptimization.mosaicSizeDegs + extraSizeDegs*[1 1];

% Add the specific targetVisualSTFdescriptor
if (isempty(targetVisualSTFdescriptor))
	rgcMosaicSurroundOptimization.targetVisualSTFdescriptor = 'default';
else
	rgcMosaicSurroundOptimization.targetVisualSTFdescriptor = targetVisualSTFdescriptor;
end


% Add all the specified surround optimization params
pStruct.rgcMosaicSurroundOptimization = rgcMosaicSurroundOptimization;


if (1==2)
	% Test to validate against the RFs depicted in Fig 4 of the Field et al (2010) paper,
	% In the supplementary section of that paper it is stated that all cells were recorded at
	% an eccentricity of 6.75 mm temporal eccentricity along the horizontal meridian.
	% According to the following computation, 6.75 mm in macaque retina corresponds to 25 degs in human retina
	Field2010Figure4CellsEccentricityMMs = mRGCMosaic.eccentricityMMsForFigure4OfFieldEtAl2010

	% Convert to degs in macaque retina
	Field2010Figure4CellsEccentricityDegs = mRGCMosaic.eccentricityDegsMacaqueRetinaForFigure4OfFieldEtAl2010

	% Convert to degs in human retina
	equivalentEccDegsInHumanRetina = mRGCMosaic.eccentricityDegsHumanRetinaForFigure4OfFieldEtAl2010
end


end