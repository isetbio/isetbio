%
% RGCMosaicAnalyzer.compute.photocurrentsForOneStimulusPeriod
%
function [temporalSupportPhotocurrent, theConePhotocurrents, theConeBackgroundPhotocurrents, ...
          thePeriodicConeMosaicExcitations, thePeriodicConeMosaicExcitationsTemporalSupport] = ...
          photocurrentsForOneStimulusPeriod(...
            eccentricityDegs, temporalSupportSeconds, theConeMosaicExcitationResponseSequence, ...
            nWarmUpPeriods, pCurrentTemporalResolutionSeconds, osTimeStep, coneTypes, ...
            debugInputConeMosaicPcurrentResponse, plotTitle, varargin)

    p = inputParser;
    p.addParameter('onlyKeepResponseDuringLastStimulusPeriod', true, @islogical);
    p.addParameter('computePhotocurrentResponsesOnlyForSelectConeIndices', [], @isnumeric);
    % Execute the parser
    p.parse(varargin{:});
    onlyKeepResponseDuringLastStimulusPeriod = p.Results.onlyKeepResponseDuringLastStimulusPeriod;
    computePhotocurrentResponsesOnlyForSelectConeIndices = p.Results.computePhotocurrentResponsesOnlyForSelectConeIndices;

    nTimeBins = size(theConeMosaicExcitationResponseSequence,1);
    nCones = size(theConeMosaicExcitationResponseSequence,2);
    dT = temporalSupportSeconds(2)-temporalSupportSeconds(1);
    periodSeconds = temporalSupportSeconds(nTimeBins)+dT;
    temporalSupportSecondsPeriodic = reshape(temporalSupportSeconds, [1 numel(temporalSupportSeconds)]);
    for k = 1:nWarmUpPeriods
        temporalSupportSecondsPeriodic = cat(2, temporalSupportSecondsPeriodic, temporalSupportSeconds+k*periodSeconds);
    end
    
    % Cone mosaic integration time
    coneMosaicIntegrationTime = temporalSupportSeconds(2)-temporalSupportSeconds(1);
    
    % Preallocate memory
    theSingleConeExcitations = theConeMosaicExcitationResponseSequence(:,1);
 

    temporalSupportPhotocurrent = computeSingleConePhotocurrentResponse( ...
        eccentricityDegs, theSingleConeExcitations, temporalSupportSeconds, ...
        nWarmUpPeriods, onlyKeepResponseDuringLastStimulusPeriod, ...
        coneMosaicIntegrationTime, osTimeStep, pCurrentTemporalResolutionSeconds, false);


    theConePhotocurrents = zeros(numel(temporalSupportPhotocurrent), nCones);
    theConeBackgroundPhotocurrents = zeros(1, nCones);
    
    if (debugInputConeMosaicPcurrentResponse)

        thePeriodicConeMosaicExcitations = [];

        % Sort cones according to responses
        [~,sortedConeIndices] = sort(max(abs(theConeMosaicExcitationResponseSequence), [],1), 'descend');

        visualizedConeIndices = [];
        allConeIndicesForDebugVisualization = sortedConeIndices(1:100);

        if (~isempty(allConeIndicesForDebugVisualization))
            % We will visualize cone excitations and photocurrents for 3 cones, 1L, 1M, 1S
            idxL = find(coneTypes(allConeIndicesForDebugVisualization) == cMosaic.LCONE_ID);
            idxM = find(coneTypes(allConeIndicesForDebugVisualization) == cMosaic.MCONE_ID);
            idxS = find(coneTypes(allConeIndicesForDebugVisualization) == cMosaic.SCONE_ID);

            if (~isempty(idxL))
                visualizedConeIndices(numel(visualizedConeIndices)+1) = allConeIndicesForDebugVisualization(idxL(1));
            end
            if (~isempty(idxM))
                visualizedConeIndices(numel(visualizedConeIndices)+1) = allConeIndicesForDebugVisualization(idxM(1));
            end
            if (~isempty(idxS))
                visualizedConeIndices(numel(visualizedConeIndices)+1) = allConeIndicesForDebugVisualization(idxS(1));
            end

        else
            fprintf('Did not find any cones with significant response. No debug visualization')
        end

        for iDebugConeIndex = 1:numel(visualizedConeIndices)

            iConeIndex = visualizedConeIndices(iDebugConeIndex);
            fprintf('Serially computing (debugInputConeMosaicPcurrentResponse is ON) photocurrent for cone %d of %d\n', iConeIndex, size(theConeMosaicExcitationResponseSequence,2));
    
            % Retrieve the cone excitation response
            theSingleConeExcitations = theConeMosaicExcitationResponseSequence(:,iConeIndex);
   

            % Compute photocurrent for this cone making the cone excitation response periodic by concatenating nWarmUpPeriods
            [temporalSupportPhotocurrent, ...
             theConePhotoCurrentDifferentialResponse, ...
             theConeBackgroundPhotoCurrent, ...
             theConeExcitationsSingleConePeriodic, ...
             photocurrentResponseTimeAxisPeriodic, ...
             thePcurrentResponsePeriodic, ...
             thePcurrentBackgroundResponseTransient] = computeSingleConePhotocurrentResponse(...
                    eccentricityDegs, theSingleConeExcitations, temporalSupportSeconds, ...
                    nWarmUpPeriods, onlyKeepResponseDuringLastStimulusPeriod, ...
                    coneMosaicIntegrationTime, osTimeStep, pCurrentTemporalResolutionSeconds, true);
   
            % Retrieve the photocurrent response
            theConePhotocurrents(:, iConeIndex) = theConePhotoCurrentDifferentialResponse;
            theConeBackgroundPhotocurrents(iConeIndex) = theConeBackgroundPhotoCurrent;
    
            if (ismember(iConeIndex, visualizedConeIndices))
                hFig(1) = figure(33); clf;
                set(hFig, 'Position', [10 10 1600 1150], 'Color', [1 1 1]);

                % Contrast cone excitations to photocurrents
                coneExcitationsTransformedIntoResponseModulations = false;
    
                RGCMosaicAnalyzer.visualize.coneExcitationsVsPhotocurrentResponse(...
                    hFig, ...
                    plotTitle, ...
                    coneTypes(iConeIndex), ...
                    temporalSupportSeconds, ...
                    theSingleConeExcitations, ...
                    coneExcitationsTransformedIntoResponseModulations, ...
                    temporalSupportPhotocurrent, ...
                    theConePhotoCurrentDifferentialResponse, ...
                    theConeBackgroundPhotoCurrent, ...
                    theConeExcitationsSingleConePeriodic, ...
                    temporalSupportSecondsPeriodic, ...
                    photocurrentResponseTimeAxisPeriodic, ...
                    thePcurrentResponsePeriodic, ...
                    thePcurrentBackgroundResponseTransient);
               
                disp('In debug mode. Hit enter to visualize another response, or Control C to exit.')
                pause();
            end
        end % for iDebugConeIndex

    else  % if (!debugInputConeMosaicPcurrentResponse)
    
        % Allocate memory
        theSingleConeExcitations = theConeMosaicExcitationResponseSequence(:,1);
        [~, ~, ~, theSingleConeExcitationPeriodicResponse] = computeSingleConePhotocurrentResponse(...
                    eccentricityDegs, theSingleConeExcitations, temporalSupportSeconds, ...
                    nWarmUpPeriods, onlyKeepResponseDuringLastStimulusPeriod, ...
                    coneMosaicIntegrationTime, osTimeStep, ...
                    pCurrentTemporalResolutionSeconds, true);
        thePeriodicConeMosaicExcitations = zeros(nCones, numel(theSingleConeExcitationPeriodicResponse));

        tic

        if (isempty(computePhotocurrentResponsesOnlyForSelectConeIndices))
            coneIndicesToComputePhotocurrentsFor = 1:nCones;
        else
            coneIndicesToComputePhotocurrentsFor = computePhotocurrentResponsesOnlyForSelectConeIndices;
        end

        parfor iConeIndex = 1:nCones

            if (mod(iConeIndex,100) == 0)
                % Give some feedback
                fprintf('Computing photocurrent for cone %d of %d\n', iConeIndex, size(theConeMosaicExcitationResponseSequence,2));
            end
    
            if (~ismember(iConeIndex, coneIndicesToComputePhotocurrentsFor))
                continue;
            end

            % Retrieve the cone excitation response
            theSingleConeExcitations = theConeMosaicExcitationResponseSequence(:,iConeIndex);

            % Compute photocurrent for this cone making the cone excitation response periodic by concatenating nWarmUpPeriods
            [~, theConePhotocurrents(:, iConeIndex), ...
                theConeBackgroundPhotocurrents(iConeIndex) , ...
                thePeriodicConeMosaicExcitations(iConeIndex,:)] = computeSingleConePhotocurrentResponse(...
                    eccentricityDegs, theSingleConeExcitations, temporalSupportSeconds, ...
                    nWarmUpPeriods, onlyKeepResponseDuringLastStimulusPeriod, ...
                    coneMosaicIntegrationTime, osTimeStep, ...
                    pCurrentTemporalResolutionSeconds, true);

        end % parfor iConeIndex
    
        fprintf('\nPhotocurrent computation took %2.1f minutes\n\n', toc/60);
    end  % if (!debugInputConeMosaicPcurrentResponse)

    thePeriodicConeMosaicExcitationsTemporalSupport = (0:size(thePeriodicConeMosaicExcitations,2)-1) * coneMosaicIntegrationTime;

    % Shape the same way as the excitations response
    thePeriodicConeMosaicExcitations = thePeriodicConeMosaicExcitations';
    thePeriodicConeMosaicExcitationsTemporalSupport = thePeriodicConeMosaicExcitationsTemporalSupport';
end


function [temporalSupportPhotocurrent, ...
    theConePhotoCurrentDifferentialResponse, ...
    backgroundPhotocurrent, ...
    theConeExcitationsPeriodic, ...
    photocurrentPeriodicResponseTimeAxis, ...
    thePcurrentDifferentialPeriodicResponse, ...
    thePcurrentBackgroundResponseTransient] = computeSingleConePhotocurrentResponse(...
        eccentricityDegs, theConeExcitations, temporalSupportSeconds, ...
        nWarmUpPeriods, onlyKeepResponseDuringLastStimulusPeriod, ...
        coneMosaicIntegrationTime, osTimeStep, pCurrentTemporalResolutionSeconds, ...
        skipAssertions)

% Make the cone excitation response periodic by concatenating nWarmUpPeriods
theConeExcitations = theConeExcitations(1:end);
theConeExcitationsPeriodic = reshape(theConeExcitations, [1 numel(theConeExcitations)]);
for k = 1:nWarmUpPeriods
    theConeExcitationsPeriodic = cat(2, theConeExcitationsPeriodic, theConeExcitationsPeriodic(1:numel(theConeExcitations)));
end


% Convert excitation counts to excitation rates
theConeExcitationsRatePeriodic = theConeExcitationsPeriodic(:) / coneMosaicIntegrationTime;
backgroundConeExcitationRate = mean(theConeExcitationsRatePeriodic);

% Compute the pCurrent response to the periodic stimulus
[thePcurrentDifferentialPeriodicResponse, photocurrentPeriodicResponseTimeAxis, thePcurrentBackgroundResponseTransient] = ...
    cMosaic.photocurrentFromConeExcitationRateUsingBiophysicalOSmodel(sqrt(sum(eccentricityDegs(:).^2)), ...
        theConeExcitationsRatePeriodic, ...
        backgroundConeExcitationRate, ...
        coneMosaicIntegrationTime, ...
        pCurrentTemporalResolutionSeconds, ...
        'osTimeStepSeconds', osTimeStep, ...
        'skipAssertions', skipAssertions);

backgroundPhotocurrent = thePcurrentBackgroundResponseTransient(end);

if (onlyKeepResponseDuringLastStimulusPeriod)
    % Only keep pCurrent response during the last stimulus period
    dToriginal = temporalSupportSeconds(2)-temporalSupportSeconds(1);
    dT = photocurrentPeriodicResponseTimeAxis(2)-photocurrentPeriodicResponseTimeAxis(1);
    tOneStimulusCycle = temporalSupportSeconds(end)-temporalSupportSeconds(1);
    idx = find(photocurrentPeriodicResponseTimeAxis >= photocurrentPeriodicResponseTimeAxis(end)-(tOneStimulusCycle+0.5*dToriginal-dT));

    theConePhotoCurrentDifferentialResponse = thePcurrentDifferentialPeriodicResponse(idx);
    temporalSupportPhotocurrent = photocurrentPeriodicResponseTimeAxis(idx);
else
    theConePhotoCurrentDifferentialResponse = thePcurrentDifferentialPeriodicResponse;
    temporalSupportPhotocurrent = photocurrentPeriodicResponseTimeAxis;
end


% Time support starts at 0 msec
temporalSupportPhotocurrent = temporalSupportPhotocurrent - temporalSupportPhotocurrent(1);
temporalSupportPhotocurrent = temporalSupportPhotocurrent';
end