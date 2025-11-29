function t_mRGCMosaicNonLinearities(options)

% History:
%    11/19/25  NPC  Wrote it.

% Examples:
%{

%
% NOTE: To run any RGC-related ISETBio code, such as this tutorial, users must follow
% the directions discribed in:
%    https://github.com/isetbio/isetbio/wiki/Retinal-ganglion-cell-(RGC)-mosaics
% under section "Configuring ISETBio to access RGC resources and run RGC simulations"
%

    % Basic call
    t_mRGCMosaicNonLinearities();

    t_mRGCMosaicNonLinearities( ...
        'visualizeStimulusSequence', true, ...
        'visualizeMosaicResponses', true);

    backgroundLuminanceCdM2 = 200;
    chromaticity = 'Achromatic';
    contrast = 0.5;
    temporalFrequencyHz = 0.2;
    orientationsDegs = 0;
    spatialFrequencyCPD = 3.0;

    t_mRGCMosaicNonLinearities( ...
        'backgroundLuminanceCdM2', backgroundLuminanceCdM2, ...
        'chromaticity', chromaticity, ...
        'contrast', contrast, ...
        'temporalFrequencyHz', temporalFrequencyHz, ...
        'orientationsDegs', orientationsDegs, ...
        'spatialFrequencyCPD', spatialFrequencyCPD, ...
        'computeInputConeMosaicResponses', true, ...
        'computeInputConeMosaicResponsesBasedOnConeExcitations', ~true, ...
        'computeInputConeMosaicResponsesBasedOnPhotocurrents', ~true, ...
        'computeMRGCMosaicResponses', ~true);


    coneExcitationsResponseBasedNonLinearitiesList{1} = struct(...
        'sourceSignal', 'compositeResponse', ...    % where to apply the non-linerity. choose from {'centerComponentResponse', 'surroundComponentResponse', 'compositeResponse'}
        'type', 'Naka Rushton', ...
        'params', struct(...
            'rectification', 'half', ...            % apply non-linearity to both response polarities % choose between {'half', 'full', 'none'}
            'n', 1.4, ...                           % exponent
            's', 1.0, ...                           % super-saturation exponent (super-saturating if s > 1)
            'bias', 0.0, ...                        % bias (in normalized range, based on the passed maxLinearResponse)
            'c50', 0.07, ...                         % semi-saturation response (in normalized range, based on the passed maxLinearResponse)
            'gain', 0.1, ...                        % post-non-linearity gain
            'maxLinearResponse', 1.0) ...           % max absolute value that the linear excitations-based response can have (clipping after that)
        );

    photocurrentsResponseBasedNonLinearitiesList{1} = struct(...
        'sourceSignal', 'compositeResponse', ...    % where to apply the non-linerity. choose from {'centerComponentResponse', 'surroundComponentResponse', 'compositeResponse'}
        'type', 'Naka Rushton', ...
        'params', struct(...
            'rectification', 'none',  ...           % apply non-linearity to both response polarities % choose between {'half', 'full', 'none'}
            'n', 1.4, ...                           % exponent
            's', 1.0, ...                           % super-saturation exponent (super-saturating if s > 1)
            'bias', 0.0, ...                        % bias (in normalized range, based on the passed maxLinearResponse)
            'c50', 0.07, ...                         % semi-saturation response (in normalized range, based on the passed maxLinearResponse)
            'gain', 0.1, ...                        % post-non-linearity gain
            'maxLinearResponse', 50.0) ...          % max absolute value that the linear photocurrent-based response can have (clipping after that)
        );


    coneExcitationsBasedNonLinearityParamsStruct = struct(...
        'label', 'noSatNL', ...     % some string encoding the non-linear processing; it is encoded in the output filename
        'nonLinearitiesList', coneExcitationsResponseBasedNonLinearitiesList ...
    );

    photocurrentsBasedNonLinearityParamsStruct = struct(...
        'label', coneExcitationsBasedNonLinearityParamsStruct.label, ...
        'nonLinearitiesList', photocurrentsResponseBasedNonLinearitiesList ...
    );


    t_mRGCMosaicNonLinearities( ...
        'backgroundLuminanceCdM2', backgroundLuminanceCdM2, ...
        'chromaticity', chromaticity, ...
        'contrast', contrast, ...
        'temporalFrequencyHz', temporalFrequencyHz, ...
        'orientationsDegs', orientationsDegs, ...
        'spatialFrequencyCPD', spatialFrequencyCPD, ...
        'computeInputConeMosaicResponses', true, ...
        'computeInputConeMosaicResponsesBasedOnConeExcitations', true, ...
        'computeInputConeMosaicResponsesBasedOnPhotocurrents', ~true, ...
        'computeMRGCMosaicResponses', ~true, ...
        'coneExcitationsBasedNonLinearityParamsStruct', coneExcitationsBasedNonLinearityParamsStruct, ...
        'photocurrentsBasedNonLinearityParamsStruct', photocurrentsBasedNonLinearityParamsStruct);


%}

arguments

    % ---- Mosaic specifiers for selecting a prebaked mRGC mosaic ------
    
    % See RGCMosaicConstructor.helper.utils.initializeRGCMosaicGenerationParameters
    % for what is available and to add new mosaics
    options.rgcMosaicName (1,:) char = 'PLOSpaperNasal2DegsTinyMosaic';


    % ---- Which species to employ ----
    % Choose between {'macaque', 'human'}. If 'macaque' is chosen, the input
    % cone mosaic has a 1:1 L/M cone ratio.
    options.coneMosaicSpecies  (1,:) char {mustBeMember(options.coneMosaicSpecies,{'human','macaque'})} = 'human';


    % ----- Which subject optics to employ -----
    options.opticsSubjectName (1,:) ...
        char ...
        {...
        mustBeMember(options.opticsSubjectName, ...
            { ...
            'PLOSpaperDefaultSubject' ...
            'PLOSpaperSecondSubject' ...
            'VSS2024TalkFirstSubject' ...
            'VSS2024TalkSecondSubject' ...
            'PLOSpaperStrehlRatio_0.87' ...
            'PLOSpaperStrehlRatio_0.72' ...
            'PLOSpaperStrehlRatio_0.59' ...
            'PLOSpaperStrehlRatio_0.60' ...
            'PLOSpaperStrehlRatio_0.27' ...
            'PLOSpaperStrehlRatio_0.23' ...
            'PLOSpaperStrehlRatio_0.21' ...
            'PLOSpaperStrehlRatio_0.19' ...
            'PLOSpaperStrehlRatio_0.09' ...
            } ...
            ) ...
        } ...
        = 'PLOSpaperDefaultSubject';


    % ------ targetVisualSTF options ----
    % Options are : {'default', 'x1.3 RsRcRatio'}
    % These are with respect to the macaque data of the Croner & Kaplan '95 study
    % 'default': target the mean Rs/Rc, and the mean Ks/Kc (Rs/Rc)^2
    % See RGCMosaicConstructor.helper.surroundPoolingOptimizerEngine.generateTargetVisualSTFmodifiersStruct
    % for all existing options
    options.targetVisualSTFdescriptor (1,:) char = 'default';

    % Submosaic to map
    options.cropParams = [];

    % Different options for the optics
    options.opticsForResponses = [];

    % Luminance, chromaticity and contrast
    options.backgroundLuminanceCdM2 (1,1) double = 100;
    options.backgroundChromaticity (1,2) double = [0.31, 0.31];
    options.chromaticity (1,:) char{mustBeMember(options.chromaticity, {'LconeIsolating' 'MconeIsolating' 'SconeIsolating', 'Achromatic'})} = 'Achromatic';
    options.contrast (1,1) double = 1.0;

    % Temporal frequency
    options.temporalFrequencyHz (1,1) double = 4;
    
    % Decreasing the spatial phase increment results in higher temporal resolution, and vice versa
    options.spatialPhaseIncrementDegs (1,1) double = 30;

    % Orientation and spatial frequency
    options.orientationsDegs (1,1) double = 0;
    options.spatialFrequencyCPD (1,1) double = 3.0;

    % Nonlinearities
    options.coneExcitationsBasedNonLinearityParamsStruct = [];
    options.photocurrentsBasedNonLinearityParamsStruct = [];

    
    % Visualizations
    options.visualizeStimulusSequence (1,1) logical = false;
    options.visualizeMosaicResponses (1,1) logical = false;
    options.visualizeConeExcitationVsPhotocurrentBasedResponses (1,1) logical = false;

    % ---- Choices of actions to perform ----
    % Whether to compute the input cone mosaic STF responses
    options.computeInputConeMosaicResponses (1,1) logical = true;
    options.computeInputConeMosaicResponsesBasedOnConeExcitations (1,1) logical = true;
    options.computeInputConeMosaicResponsesBasedOnPhotocurrents (1,1) logical = true;

    % Whether to compute the input cone mosaic STF responses
    options.computeMRGCMosaicResponses (1,1) logical = false;
    
    % Whether to close previously open figures
    options.closePreviouslyOpenFigures (1,1) logical = true;

end % arguments



% Set flags from key/value pairs

% Mosaic specifiers for selecting a previously synthesized center-connected mRGC mosaic 
rgcMosaicName = options.rgcMosaicName;
coneMosaicSpecies = options.coneMosaicSpecies;
opticsSubjectName = options.opticsSubjectName;
targetVisualSTFdescriptor = options.targetVisualSTFdescriptor;

% Optics to employ for the computations
opticsForResponses = options.opticsForResponses;

% Luminance, chromaticity and contrast
backgroundLuminanceCdM2 = options.backgroundLuminanceCdM2;
backgroundChromaticity = options.backgroundChromaticity;
chromaticity = options.chromaticity; 
contrast = options.contrast;

% Temporal frequency
temporalFrequencyHz = options.temporalFrequencyHz;

% Spatial params
spatialPhaseIncrementDegs = options.spatialPhaseIncrementDegs;
orientationsDegs = options.orientationsDegs;
spatialFrequencyCPD = options.spatialFrequencyCPD;

% Nonlinearities
coneExcitationsBasedNonLinearityParamsStruct = options.coneExcitationsBasedNonLinearityParamsStruct;
photocurrentsBasedNonLinearityParamsStruct = options.photocurrentsBasedNonLinearityParamsStruct;


% Mosaic cropping
cropParams = options.cropParams;

% Visualizations
visualizeStimulusSequence = options.visualizeStimulusSequence;
visualizeMosaicResponses = options.visualizeMosaicResponses;
visualizeConeExcitationVsPhotocurrentBasedResponses = options.visualizeConeExcitationVsPhotocurrentBasedResponses;

% Actions to perform
computeInputConeMosaicResponses = options.computeInputConeMosaicResponses;
computeInputConeMosaicResponsesBasedOnConeExcitations = options.computeInputConeMosaicResponsesBasedOnConeExcitations;
computeInputConeMosaicResponsesBasedOnPhotocurrents = options.computeInputConeMosaicResponsesBasedOnPhotocurrents;
computeMRGCMosaicResponses = options.computeMRGCMosaicResponses;


% Load the mRGCmosaic specified by the passed parameters:
% coneMosaicSpecies, opticsSubjectName, rgcMosaicName, targetVisualSTFdescriptor
% and generate the optics that were used to synthesize the mosaic
[theMRGCmosaic, theOI, thePSFatTheMosaicEccentricity, ...
    prebakedMRGCMosaicDir, prebakedMRGCMosaicFilename] = mRGCMosaic.loadPrebakedMosaic(...
        coneMosaicSpecies, opticsSubjectName, rgcMosaicName, targetVisualSTFdescriptor, ...
        'computeTheMosaicOptics', true, ...
        'opticsToEmploy', opticsForResponses, ...
        'cropParams', cropParams);


% Filenames for intermediate responses
intermediateDataDir = RGCMosaicConstructor.filepathFor.intermediateDataDir();

% Directory of intermediateDataDir where all computed data will go
subDir = 'nonLinearityDemos';

figureDir = fullfile(isetbioRootPath,'local',mfilename);
if (~exist(figureDir,'dir'))
    mkdir(figureDir);
end
fprintf('Will save figures/videos into %s\n',figureDir);


[theInputConeMosaicResponsesFullFileName, theMRGCMosaicResponsesFullFileName] = computeSimulationForCurrentStimulus(...
    orientationsDegs, spatialFrequencyCPD, ...
    chromaticity, contrast, ...
    backgroundChromaticity, backgroundLuminanceCdM2, ...
    spatialPhaseIncrementDegs, temporalFrequencyHz, ...
    prebakedMRGCMosaicFilename, intermediateDataDir, subDir, figureDir, ...
    opticsForResponses, theOI, thePSFatTheMosaicEccentricity, ...
    theMRGCmosaic, ...
    coneExcitationsBasedNonLinearityParamsStruct, ...
    photocurrentsBasedNonLinearityParamsStruct, ...
    computeInputConeMosaicResponses, ...
    computeInputConeMosaicResponsesBasedOnConeExcitations, ...
    computeInputConeMosaicResponsesBasedOnPhotocurrents, ...
    computeMRGCMosaicResponses, ...
    visualizeMosaicResponses,visualizeStimulusSequence);


end


function [theInputConeMosaicResponsesFullFileName, theMRGCMosaicResponsesFullFileName] = computeSimulationForCurrentStimulus(...
    orientationsDegs, spatialFrequencyCPD, ...
    chromaticity, contrast, ...
    backgroundChromaticity, backgroundLuminanceCdM2, ...
    spatialPhaseIncrementDegs, temporalFrequencyHz, ...
    prebakedMRGCMosaicFilename, intermediateDataDir, subDir, figureDir, ...
    opticsForResponses, theOI, thePSFatTheMosaicEccentricity, ...
    theMRGCmosaic, ...
    coneExcitationsBasedMRGCNonLinearityParamsStruct, ...
    photocurrentsBasedMRGCNonLinearityParamsStruct, ...
    computeInputConeMosaicResponses, ...
    computeInputConeMosaicResponsesBasedOnConeExcitations, ...
    computeInputConeMosaicResponsesBasedOnPhotocurrents, ...
    computeMRGCMosaicResponses, ...
    visualizeMosaicResponses,visualizeStimulusSequence)

    % Encode chromaticity, eccentricity and size of the mosaic
    postFix = sprintf('%s@%2.0f_%2.0fCd_%.0fDeg_%.1fCpD_%2.0fHz_%s_Ecc_%2.1f_%2.1f_Size_%2.1f_%2.1f', ...
        chromaticity, ...
        contrast*100, ...
        backgroundLuminanceCdM2, ...
        orientationsDegs, ...
        spatialFrequencyCPD, ...
        temporalFrequencyHz, ...
        strrep(prebakedMRGCMosaicFilename, '.mat', ''), ...
        theMRGCmosaic.eccentricityDegs(1), ...
        theMRGCmosaic.eccentricityDegs(2), ...
        theMRGCmosaic.sizeDegs(1), ...
        theMRGCmosaic.sizeDegs(2));
   

    theInputConeMosaicResponsesFileName = fullfile(subDir, sprintf('coneResponses_%s.mat', postFix));
    theInputConeMosaicResponsesFullFileName = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
            intermediateDataDir, theInputConeMosaicResponsesFileName, ...
            'generateMissingSubDirs', true);

    theMRGCMosaicResponsesFileName = fullfile(subDir, sprintf('mRGCResponses_%s.mat', postFix));
    theMRGCMosaicResponsesFullFileName = RGCMosaicConstructor.filepathFor.augmentedPathWithSubdirs(...
            intermediateDataDir, theMRGCMosaicResponsesFileName, ...
            'generateMissingSubDirs', true);


    if (~isempty(coneExcitationsBasedMRGCNonLinearityParamsStruct))
        theMRGCMosaicResponsesFullFileName = ...
            strrep(theMRGCMosaicResponsesFullFileName, 'Responses', sprintf('Responses_%s', coneExcitationsBasedMRGCNonLinearityParamsStruct.label));
    elseif (~isempty(photocurrentsBasedMRGCNonLinearityParamsStruct))
        theMRGCMosaicResponsesFullFileName = ...
            strrep(theMRGCMosaicResponsesFullFileName, 'Responses', sprintf('Responses_%s', photocurrentsBasedMRGCNonLinearityParamsStruct.label));
    end


    % Compute input cone mosaic response
    if (computeInputConeMosaicResponses) && (computeInputConeMosaicResponsesBasedOnConeExcitations)

        % Form stimParams struct
        switch (chromaticity)
            case 'LconeIsolating'
                coneContrastDirection = [1 0 0];
            case 'MconeIsolating'
                coneContrastDirection = [0 1 0];
            case 'SconeIsolating'
                coneContrastDirection = [0 0 1];
            case 'Achromatic'
                coneContrastDirection = [1 1 1];
        end

        opticsIsDiffractionLimited = false;
        if (~isempty(opticsForResponses)) && (isstruct(opticsForResponses) && (isfield(opticsForResponses, 'type'))) && ...
           ((strcmp(opticsForResponses.type, 'adaptiveOptics6MM')) || (strcmp(opticsForResponses.type, 'adaptiveOptics6MMwithLCA')))
                opticsIsDiffractionLimited = true;
        end

        stimParams = struct(...
            'displayType', 'CRT-Sony-HorwitzLab', ...
            'displayLuminanceHeadroomPercentage', 0.5, ...
            'backgroundChromaticity', backgroundChromaticity, ...
            'backgroundLuminanceCdM2', backgroundLuminanceCdM2, ...
            'contrast', contrast, ...
            'coneContrasts', coneContrastDirection, ...
            'sizeDegs', 1.1*max(theMRGCmosaic.inputConeMosaic.sizeDegs), ...
            'spatialPhaseIncrementDegs', spatialPhaseIncrementDegs, ...
            'temporalFrequencyHz', temporalFrequencyHz, ...
            'durationSeconds', 2.0/temporalFrequencyHz, ...
            'orientationDegs', orientationsDegs, ...
            'spatialFrequencyCPD', spatialFrequencyCPD, ...
            'opticsIsDiffractionLimited', opticsIsDiffractionLimited);

        % Compute the inputConeMosaic excitation responses
        computeInputConeMosaicResponse(theMRGCmosaic, theOI, thePSFatTheMosaicEccentricity, stimParams, ...
            visualizeMosaicResponses,visualizeStimulusSequence, ...
            theInputConeMosaicResponsesFullFileName);
    end

    if (computeInputConeMosaicResponses) && (computeInputConeMosaicResponsesBasedOnPhotocurrents)

        photocurrentParams.osBiophysicalModelWarmUpTimeSeconds = 1.0;
        photocurrentParams.osBiophysicalModelTemporalResolutionSeconds = 1e-5;
        photocurrentParams.temporalResolutionSeconds =  5/1000;

        computeInputConeMosaicPhotocurrentResponse(theInputConeMosaicResponsesFullFileName, photocurrentParams);
    end


    if (computeMRGCMosaicResponses)

        if (~isempty(coneExcitationsBasedMRGCNonLinearityParamsStruct))
            coneExcitationsResponseBasedNonLinearitiesList = coneExcitationsBasedMRGCNonLinearityParamsStruct.nonLinearitiesList;
        else
            coneExcitationsResponseBasedNonLinearitiesList = {};
        end

        if (~isempty(photocurrentsBasedMRGCNonLinearityParamsStruct))
            photocurrentResponseBasedNonLinearitiesList = photocurrentsBasedMRGCNonLinearityParamsStruct.nonLinearitiesList;
        else
            photocurrentResponseBasedNonLinearitiesList = {};
        end

        computeAllMRGCMosaicResponses(...
            theInputConeMosaicResponsesFullFileName, ...
            theMRGCMosaicResponsesFullFileName, ...
            coneExcitationsResponseBasedNonLinearitiesList, ...
            photocurrentResponseBasedNonLinearitiesList, ...
            figureDir);
    end
end % computeSimulationForCurrentStimulus


%
% Supporting functions
%
%
function computeAllMRGCMosaicResponses(theInputConeMosaicResponsesFullFileName, theMRGCMosaicResponsesFullFileName, ...
    coneExcitationsResponseBasedNonLinearitiesList, ...
    photocurrentResponseBasedNonLinearitiesList, ...
    figureDir)

    % The sincle period cone excitations response of the input cone mosaic
    load(theInputConeMosaicResponsesFullFileName, ...
        'theMRGCmosaic', ...
        'stimParams', ...
        'theInputConeMosaicExcitationsResponse', ...
        'theInputConeMosaicExcitationsNullResponse');

    % Compute the NONLINEAR cone-excitations based mRGC mosaic responses to a single period
    [theMRGCmosaicExcitationBasedResponses, temporalSupportSeconds, theMRGCmosaicExcitationBasedLinearResponses] = ...
        computeMRGCMosaicResponses(theMRGCmosaic, ...
            stimParams.temporalSupportSeconds, ...
            theInputConeMosaicExcitationsResponse, ...
            coneExcitationsResponseBasedNonLinearitiesList);

    fprintf('Max LINEAR response: %f\n', max(abs(theMRGCmosaicExcitationBasedLinearResponses(:))));
    fprintf('Max NON-LINEAR response: %f\n', max(abs(theMRGCmosaicExcitationBasedResponses(:))));


    theMRGCmosaicExcitationBasedResponsesSinglePeriod = struct(...
        'linearResponses', theMRGCmosaicExcitationBasedLinearResponses, ...
        'responses', theMRGCmosaicExcitationBasedResponses, ...
        'temporalSupportSeconds', temporalSupportSeconds);

    % Save the single-period, excitations-based responses
    save(theMRGCMosaicResponsesFullFileName, ...
        'theMRGCmosaicExcitationBasedResponsesSinglePeriod', ...
        '-v7.3');

    % See if we have input cone mosaic photocurrent responses and if we do
    % compute the photocurrents based mRGC mosaic responses
    allData = who('-file', theInputConeMosaicResponsesFullFileName);

    if (ismember('theInputConeMosaicPhotocurrentsResponse', allData))

        % We have photocurrent based inputConeMosaic responses,
        % so compute mRGC responses based on those as well
         load(theInputConeMosaicResponsesFullFileName, ...
            'theInputConeMosaicPeriodicExcitationsResponse', ...
            'theInputConeMosaicPeriodicExcitationsTemporalSupportSeconds', ...
            'theInputConeMosaicPeriodicPhotocurrentTemporalSupportSeconds', ...
            'theInputConeMosaicPeriodicPhotocurrentsResponse', ...
            'theInputConeMosaicBackgroundPhotocurrents');

        % excitations-based mRGC responses (periodic)
        [theMRGCmosaicResponses, temporalSupportSeconds, theMRGCmosaicLinearResponses] = ...
            computeMRGCMosaicResponses(theMRGCmosaic, ...
                theInputConeMosaicPeriodicExcitationsTemporalSupportSeconds, ...
                theInputConeMosaicPeriodicExcitationsResponse, ...
                coneExcitationsResponseBasedNonLinearitiesList);

        theMRGCmosaicExcitationsBasedResponsesPeriodic = struct(...
            'linearResponses', theMRGCmosaicLinearResponses, ...
            'responses', theMRGCmosaicResponses, ...
            'temporalSupportSeconds', temporalSupportSeconds);

        % photocurrent-based mRGC responses (periodic)
        [theMRGCmosaicResponses, temporalSupportSeconds, theMRGCmosaicLinearResponses] = ...
            computeMRGCMosaicResponses(theMRGCmosaic, ...
                theInputConeMosaicPeriodicPhotocurrentTemporalSupportSeconds, ...
                theInputConeMosaicPeriodicPhotocurrentsResponse, ...
                photocurrentResponseBasedNonLinearitiesList);

        theMRGCmosaicPhotocurrentsBasedResponsesPeriodic = struct(...
            'linearResponses', theMRGCmosaicLinearResponses, ...
            'responses', theMRGCmosaicResponses, ...
            'temporalSupportSeconds', temporalSupportSeconds);

        % Extract single period responses
        dToriginal = theMRGCmosaicExcitationBasedResponsesSinglePeriod.temporalSupportSeconds(2)-theMRGCmosaicExcitationBasedResponsesSinglePeriod.temporalSupportSeconds(1);
        dT = temporalSupportSeconds(2)-temporalSupportSeconds(1);
        tOneStimulusCycle = theMRGCmosaicExcitationBasedResponsesSinglePeriod.temporalSupportSeconds(end)-theMRGCmosaicExcitationBasedResponsesSinglePeriod.temporalSupportSeconds(1);
        idx = find(temporalSupportSeconds >= temporalSupportSeconds(end)-(tOneStimulusCycle+0.5*dToriginal-dT));

        temporalSupportSeconds = temporalSupportSeconds(idx);
        theMRGCmosaicResponses = theMRGCmosaicResponses(:, idx, :);
        theMRGCmosaicLinearResponses = theMRGCmosaicLinearResponses(:, idx, :);

        theMRGCmosaicPhotocurrentsBasedResponsesSinglePeriod = struct(...
            'linearResponses', theMRGCmosaicLinearResponses, ...
            'responses', theMRGCmosaicResponses, ...
            'temporalSupportSeconds', temporalSupportSeconds);

        % Append the periodic, excitations, and photocurrents-based responses
        save(theMRGCMosaicResponsesFullFileName, ...
            'theMRGCmosaicExcitationsBasedResponsesPeriodic', ...
            'theMRGCmosaicPhotocurrentsBasedResponsesPeriodic', ...
            'theMRGCmosaicPhotocurrentsBasedResponsesSinglePeriod', ...
            '-append');
    end


    % Visualize responses
    mRGCsNum = size(theMRGCmosaicExcitationBasedResponsesSinglePeriod.responses,3);
    maxExcitationsBasedResponse = max(abs(theMRGCmosaicExcitationBasedResponsesSinglePeriod.responses(:)));
    if (ismember('theInputConeMosaicPhotocurrentsResponse', allData))
        idx = find(theMRGCmosaicPhotocurrentsBasedResponsesPeriodic.temporalSupportSeconds>0.5);
        maxPhotocurrentsBasedResponse = max(max(abs(squeeze(theMRGCmosaicPhotocurrentsBasedResponsesPeriodic.responses(1,idx,:)))));
    end

    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1984 768], 'Color', [1 1 1]);

    videoOBJ = VideoWriter(fullfile(figureDir,'mRGCnonLinearity'), 'MPEG-4');
    videoOBJ.FrameRate = 10;
    videoOBJ.Quality = 100;
    videoOBJ.open();


    excitationsBasedTicks = -1:0.2:1;
    photocurrentsBasedTicks = -20:2:20;

    for exemplarRGCindex = 1:mRGCsNum

        set(hFig, 'Name', sprintf('mRGC #%d of %d', exemplarRGCindex, mRGCsNum));

        % ----- The non-linear activation functions ----
        % Cone-excitations based response activation function
        ax = subplot('Position', [0.02 0.57 0.18 0.37]);
        theLinearResponse = theMRGCmosaicExcitationBasedResponsesSinglePeriod.linearResponses(1,:,exemplarRGCindex);
        theNonLinearResponse = theMRGCmosaicExcitationBasedResponsesSinglePeriod.responses(1,:,exemplarRGCindex);
        plotNonLinearity(ax, theLinearResponse, theNonLinearResponse, maxExcitationsBasedResponse, ...
            excitationsBasedTicks, false, sprintf('cone excitations - based\nactivation function'));

        if (ismember('theInputConeMosaicPhotocurrentsResponse', allData))
            % Photocurrents based response activation function (single period)
            ax = subplot('Position', [0.02 0.09 0.18 0.37]);
            theLinearResponse = theMRGCmosaicPhotocurrentsBasedResponsesSinglePeriod.linearResponses(1,:,exemplarRGCindex);
            theNonLinearResponse = theMRGCmosaicPhotocurrentsBasedResponsesSinglePeriod.responses(1,:,exemplarRGCindex);
            plotNonLinearity(ax, theLinearResponse, theNonLinearResponse, maxPhotocurrentsBasedResponse, ...
                photocurrentsBasedTicks, true, sprintf('photocurrents - based\nactivation function'));
        end


        % The full response
        if (ismember('theInputConeMosaicPhotocurrentsResponse', allData))
            % The periodic cone excitations based mRGC response
            ax = subplot('Position', [0.25 0.57 0.55 0.37]);
            theResponse = theMRGCmosaicExcitationsBasedResponsesPeriodic.responses(1,:,exemplarRGCindex);
            theLinearResponse = theMRGCmosaicExcitationsBasedResponsesPeriodic.linearResponses(1,:,exemplarRGCindex);
            theTemporalSupport = theMRGCmosaicExcitationsBasedResponsesPeriodic.temporalSupportSeconds;
            plotResponse(ax, theTemporalSupport, theResponse, theLinearResponse, maxExcitationsBasedResponse, ...
                excitationsBasedTicks, false, true, 'cone excitations - based mRGC responses (full simulation time)');

            % The periodic photocurrents based mRGC response
            ax = subplot('Position', [0.25 0.09 0.55 0.37]);
            theResponse = theMRGCmosaicPhotocurrentsBasedResponsesPeriodic.responses(1,:,exemplarRGCindex);
            theLinearResponse = theMRGCmosaicPhotocurrentsBasedResponsesPeriodic.linearResponses(1,:,exemplarRGCindex);
            theTemporalSupport = theMRGCmosaicPhotocurrentsBasedResponsesPeriodic.temporalSupportSeconds;
            plotResponse(ax, theTemporalSupport, theResponse, theLinearResponse, maxPhotocurrentsBasedResponse, ...
                photocurrentsBasedTicks, true, true, 'photocurrents - based mRGC responses (full simulation time)');
        end



        % One period of the cone excitations based mRGC response
        ax = subplot('Position', [0.86 0.57 0.13 0.35]);
        theLinearResponse = theMRGCmosaicExcitationBasedResponsesSinglePeriod.linearResponses(1,:,exemplarRGCindex);
        theResponse = theMRGCmosaicExcitationBasedResponsesSinglePeriod.responses(1,:,exemplarRGCindex);
        theTemporalSupport = theMRGCmosaicExcitationBasedResponsesSinglePeriod.temporalSupportSeconds;
        plotResponse(ax, theTemporalSupport, theResponse, theLinearResponse, maxExcitationsBasedResponse, ...
            excitationsBasedTicks, false, false, sprintf('1 period'));

        if (ismember('theInputConeMosaicPhotocurrentsResponse', allData))
            % One period of the photocurrents based mRGC response
            ax = subplot('Position', [0.86 0.09 0.13 0.35]);
            theResponse = theMRGCmosaicPhotocurrentsBasedResponsesSinglePeriod.responses(1,:,exemplarRGCindex);
            theLinearResponse = theMRGCmosaicPhotocurrentsBasedResponsesSinglePeriod.linearResponses(1,:,exemplarRGCindex);
            theTemporalSupport = theMRGCmosaicPhotocurrentsBasedResponsesSinglePeriod.temporalSupportSeconds;
            plotResponse(ax, theTemporalSupport, theResponse, theLinearResponse, maxPhotocurrentsBasedResponse, ...
                photocurrentsBasedTicks, true, false, sprintf('1 period'));
        end




        drawnow;
        videoOBJ.writeVideo(getframe(hFig));
    end

    videoOBJ.close();

end


function plotResponse(ax, theTemporalSupport, theResponse, theLinearResponse, maxResponse, yTicks, showXlabel, showLegend, plotTitle)
    linearColor = [1 0 0];
    nonLinearColor = [0.5 0.5 1.0];
    plot(ax, theTemporalSupport*1e3, theLinearResponse, 'k-', 'LineWidth', 2.0, 'Color', 0.75*linearColor);
    hold(ax, 'on')
    p1 = scatter(ax, theTemporalSupport*1e3, theLinearResponse, 100, ...
        'MarkerFaceColor', linearColor, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeColor', 0.75*linearColor);
    plot(ax, theTemporalSupport*1e3, theResponse, '-', 'LineWidth', 2.0, 'Color', nonLinearColor);
    p2 = scatter(ax, theTemporalSupport*1e3, theResponse, 100,  ...
        'MarkerFaceColor', nonLinearColor, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeColor', 0.75*nonLinearColor);
    hold(ax, 'on')
    plot(ax, theTemporalSupport*1e3, theResponse*0, 'k-', 'LineWidth', 1.0);
    if (showLegend)
        legend(ax, [p1 p2], {'linear', 'non-linear'}, 'Location', 'NorthEast', 'NumColumns', 1);
    end
    hold(ax, 'off')
    set(ax, 'YLim', maxResponse * [-1.05 1.05], 'XLim', [theTemporalSupport(1) theTemporalSupport(end)]*1e3, 'FontSize', 16);
    set(ax, 'YTick', yTicks, 'XTick', 0:50:10000);
    xtickangle(ax, 90);
    grid(ax, 'on')
    title(plotTitle);
    if (showXlabel)
        xlabel(ax, 'time (seconds)')
    end
    ylabel(ax, 'mRGC response');
end


function plotNonLinearity(ax, theLinearResponse, theResponse, maxResponse, xyTicks, showXlabel, plotTitle)
    plot(ax, maxResponse*[-1 1], [0 0], 'k-', 'LineWidth', 1.0);
    hold(ax, 'on')
    plot(ax, [0 0], maxResponse*[-1 1], 'k-', 'LineWidth', 1.0);
    %plot(ax, theLinearResponse, theResponse, 'r-', 'LineWidth', 1.5);

    scatter(ax, theLinearResponse, theResponse, 100, ...
        'LineWidth', 1.0, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5, ...
        'MarkerFaceColor', [0.5 0.5 1], 'MarkerEdgeColor', 'none');
    hold(ax, 'off')

    set(ax, 'YLim', maxResponse * [-1 1], 'XLim', maxResponse * [-1 1], 'FontSize', 16);
    set(ax, 'XTick', xyTicks, 'YTick', xyTicks);
    xtickangle(ax, 90);
    grid(ax, 'on')
    axis(ax, 'square');
    title(plotTitle);
    if (showXlabel)
        xlabel(ax, 'linear mRGC response');
    end
    ylabel(ax, 'mRGC response');
end


function [theMRGCmosaicResponses, theMRGCmosaicResponseTemporalSupportSeconds, theLinearMRGCmosaicResponses] = ...
    computeMRGCMosaicResponses(...
        theMRGCmosaic, ...
        theInputConeMosaicTemporalSupportSeconds, ...
        theInputConeMosaicResponse, ...
        nonLinearitiesList)

    nTimeBins = size(theInputConeMosaicResponse,1);
    nCones = size(theInputConeMosaicResponse,2);

    % Compute the mRGCmosaic response
    [theMRGCmosaicResponses, ~, theMRGCmosaicResponseTemporalSupportSeconds, theLinearMRGCmosaicResponses] = theMRGCmosaic.compute(...
            reshape(theInputConeMosaicResponse, [1 nTimeBins nCones]), ...
            theInputConeMosaicTemporalSupportSeconds, ...
            'nonLinearitiesList', nonLinearitiesList);

end


function computeInputConeMosaicPhotocurrentResponse(theInputConeMosaicResponsesFullFileName, photocurrentParams)

    load(theInputConeMosaicResponsesFullFileName, ...
        'theMRGCmosaic', ...
        'stimParams', ...
        'theInputConeMosaicExcitationsTemporalSupportSeconds', ...
        'theInputConeMosaicExcitationsResponse', ...
        'theInputConeMosaicExcitationsNullResponse');

    % Transform cone modulations to cone excitations
    if (stimParams.coneMosaicModulationBasedResponse)
        theInputConeMosaicExcitationsNullResponse = reshape(theInputConeMosaicExcitationsNullResponse, [1 numel(theInputConeMosaicExcitationsNullResponse)]);

        % Transform from modulations to excitations
        %Rmod = (Recx-Ro)/Ro -> (Rmod+1)*Ro = RExc
        theInputConeMosaicExcitationsResponse = theInputConeMosaicExcitationsResponse + 1;
        theInputConeMosaicExcitationsResponse = bsxfun(@times, theInputConeMosaicExcitationsResponse, theInputConeMosaicExcitationsNullResponse);
    end

    % Compute # of warm up periods
    stimulusPeriodDuration = 1/stimParams.temporalFrequencyHz;

    nWarmUpPeriods = ceil(photocurrentParams.osBiophysicalModelWarmUpTimeSeconds/stimulusPeriodDuration);

    debugInputConeMosaicPcurrentResponse = false;
    plotTitle = 'debug';

    % Compute photocurrents
    [theInputConeMosaicPeriodicPhotocurrentTemporalSupportSeconds, ...
     theInputConeMosaicPeriodicPhotocurrentsResponse, theInputConeMosaicBackgroundPhotocurrents, ...
     theInputConeMosaicPeriodicExcitationsResponse, ...
     theInputConeMosaicPeriodicExcitationsTemporalSupportSeconds] = RGCMosaicAnalyzer.compute.photocurrentsForOneStimulusPeriod(...
        theMRGCmosaic.eccentricityDegs, ...
        stimParams.temporalSupportSeconds(1:end-1), ...   % Get rid of last point which is a repeat of the first point
        theInputConeMosaicExcitationsResponse(1:end-1,:), ...       % Get rid of last point which is a repeat of the first point
        nWarmUpPeriods, ...
        photocurrentParams.temporalResolutionSeconds, ...
        photocurrentParams.osBiophysicalModelTemporalResolutionSeconds, ...
        theMRGCmosaic.inputConeMosaic.coneTypes, ...
        debugInputConeMosaicPcurrentResponse, ...
        plotTitle, ...
        'onlyKeepResponseDuringLastStimulusPeriod', false);

     % Transform cone excitations back to cone modulations
    if (stimParams.coneMosaicModulationBasedResponse)
        % Transform from modulations to excitations
        %Rmod = (Recx-Ro)/Ro -> (Rmod+1)*Ro = RExc

        idx = find(theInputConeMosaicExcitationsNullResponse == 0);
        normalizingResponse = 1./ theInputConeMosaicExcitationsNullResponse;
        normalizingResponse(idx) = 0;
        theInputConeMosaicPeriodicExcitationsResponse = bsxfun(@times, theInputConeMosaicPeriodicExcitationsResponse, normalizingResponse) - 1;
    end


    save(theInputConeMosaicResponsesFullFileName, ...
        'theInputConeMosaicPeriodicPhotocurrentTemporalSupportSeconds', ...
        'theInputConeMosaicPeriodicPhotocurrentsResponse', ...
        'theInputConeMosaicBackgroundPhotocurrents', ...
        'theInputConeMosaicPeriodicExcitationsResponse', ...
        'theInputConeMosaicPeriodicExcitationsTemporalSupportSeconds', ...
        '-append');
end


function computeInputConeMosaicResponse(theMRGCmosaic, theOI, thePSFatTheMosaicEccentricity, stimParams, ...
    visualizeResponse, visualizeStimulusSequence, ...
    theInputConeMosaicResponsesFullFileName)

    % Determine the stimulus pixel resolution to be a fraction of the minimum cone aperture 
    % or cone spacing in the mosaic here, half of the cone spacing
    theMetric = 'cone aperture';  % choose from {'cone aperture' or cone spacing'}
    if (stimParams.opticsIsDiffractionLimited)
        theFraction = 0.1;
    else
        theFraction = 0.25;
    end
    targetRGCindices =  1:theMRGCmosaic.rgcsNum;
    stimulusResolutionDegs = ...
        RGCMosaicConstructor.helper.simulateExperiment.stimulusResolutionFromConeApertureOrConeSpacing(...
            theMRGCmosaic, targetRGCindices, theFraction, theMetric);
    stimParams.resolutionDegs = stimulusResolutionDegs;

    % Generate presentation display, 20% luminance headroom
    viewingDistanceMeters = 4;
    
    % Generate presentation display
    thePresentationDisplay = visualStimulusGenerator.presentationDisplay(...
        theMRGCmosaic.inputConeMosaic.wave, ...
        stimulusResolutionDegs, ...
        viewingDistanceMeters, ...
        'displayType', stimParams.displayType, ...
        'meanLuminanceCdPerM2', stimParams.backgroundLuminanceCdM2, ...
        'luminanceHeadroom', stimParams.displayLuminanceHeadroomPercentage, ...
        'bitDepth', 20);


    % Generate the spatial modulation patterns for all spatial phases of the drifting grating
    [theDriftingGratingSpatialModulationPatterns, ...
     spatialSupportDegs, spatialPhasesDegs, ...
     temporalSupportSeconds, temporalRamp] = visualStimulusGenerator.driftingGratingModulationPatterns(stimParams);


    % Generate scenes for the different frames of the drifting grating and for the null stimulus
    [theDriftingGratingFrameScenes, theNullStimulusScene] = visualStimulusGenerator.stimulusFramesScenes(...
        thePresentationDisplay, stimParams, theDriftingGratingSpatialModulationPatterns, ...
        'frameIndexToCompute', [], ... % [] field indicates that all stimulus frame scenes must be computed
        'validateScenes', ~true);

    % Compute input cone mosaic response to this orientation & spatial frequency
    stimulusPosition = 'mosaic-centered';

    % Set the mosaic integration time equal to the duration of one stimulus frame
    theMRGCmosaic.inputConeMosaic.integrationTime = temporalSupportSeconds(2)-temporalSupportSeconds(1);

    % Cone-modulation based responses
    stimParams.coneMosaicModulationBasedResponse = true;

    % Compute cone mosaic responses to each stimulus frame
    [theInputConeMosaicExcitationsResponse, theInputConeMosaicExcitationsNullResponse] = ...
        RGCMosaicConstructor.helper.simulateExperiment.inputConeMosaicResponseToStimulusFrameSequence(...
            theMRGCmosaic, ...
            theOI, ...
            theNullStimulusScene, ...
            theDriftingGratingFrameScenes, ...
            stimulusPosition, ...
            stimParams.coneMosaicModulationBasedResponse, ...
            'visualizeResponse', visualizeResponse, ...
            'visualizeStimulusSequence', visualizeStimulusSequence, ...
            'thePresentationDisplayForVisualizingOpticalSceneOrImage', thePresentationDisplay);

    % Temporal support for cone mosaic excitations
    theInputConeMosaicExcitationsTemporalSupportSeconds = temporalSupportSeconds;

    % Update stimParams
    stimParams.spatialPhasesDegs = spatialPhasesDegs;
    stimParams.temporalSupportSeconds = temporalSupportSeconds;
    stimParams.temporalRamp = temporalRamp;

    % Save results
    fprintf('Saving computed input cone mosaic STF cone excitation responses to %s.\n', theInputConeMosaicResponsesFullFileName);
    save(theInputConeMosaicResponsesFullFileName, ...
        'theMRGCmosaic', ...
        'thePSFatTheMosaicEccentricity', ...
        'stimParams', ...
        'theInputConeMosaicExcitationsTemporalSupportSeconds', ...
        'theInputConeMosaicExcitationsResponse', ...
        'theInputConeMosaicExcitationsNullResponse', ...
        '-v7.3');

end

