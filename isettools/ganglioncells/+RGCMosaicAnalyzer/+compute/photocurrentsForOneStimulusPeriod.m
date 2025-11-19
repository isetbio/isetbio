%
% RGCMosaicAnalyzer.compute.photocurrentsForOneStimulusPeriod
%
function [temporalSupportPhotocurrent, theConePhotocurrents, theConeBackgroundPhotocurrents] = ...
    photocurrentsForOneStimulusPeriod(eccentricityDegs, temporalSupportSeconds, theConeMosaicExcitationResponseSequence, ...
    nWarmUpPeriods, pCurrentTemporalResolutionSeconds, osTimeStep, coneTypes, debugInputConeMosaicPcurrentResponse, plotTitle)

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
        eccentricityDegs, theSingleConeExcitations, temporalSupportSeconds, nWarmUpPeriods, ...
        coneMosaicIntegrationTime, osTimeStep, pCurrentTemporalResolutionSeconds, false);
    
    theConePhotocurrents = zeros(numel(temporalSupportPhotocurrent), nCones);
    theConeBackgroundPhotocurrents = zeros(1, nCones);
    
    if (debugInputConeMosaicPcurrentResponse)
        % We will visualize cone excitations and photocurrents for 3 cones, 1L, 1M, 1S
        idxL = find(coneTypes == cMosaic.LCONE_ID);
        idxM = find(coneTypes == cMosaic.MCONE_ID);
        idxS = find(coneTypes == cMosaic.SCONE_ID);
    
        visualizedConeIndices = [];
        if (~isempty(idxL))
            visualizedConeIndices(numel(visualizedConeIndices)+1) = idxL(1);
        end
        if (~isempty(idxM))
            visualizedConeIndices(numel(visualizedConeIndices)+1) = idxM(1);
        end
        if (~isempty(idxS))
            visualizedConeIndices(numel(visualizedConeIndices)+1) = idxS(1);
        end
    
        hFig(1) = figure(33); clf;
        set(hFig, 'Position', [10 10 1600 1150], 'Color', [1 1 1]);

        for iConeIndex = 1:size(theConeMosaicExcitationResponseSequence,2)
            fprintf('Computing photocurrent for cone %d of %d\n', iConeIndex, size(theConeMosaicExcitationResponseSequence,2));
    
            % Retrieve the cone excitation response
            theSingleConeExcitations = theConeMosaicExcitationResponseSequence(:,iConeIndex);
    
            % Compute photocurrent for this cone making the cone excitation response periodic by concatenating nWarmUpPeriods
            [~, theConePhotoCurrentDifferentialResponse, ...
                theConeBackgroundPhotoCurrent, ...
                theConeExcitationsSingleConePeriodic, ...
                photocurrentResponseTimeAxisPeriodic, ...
                thePcurrentResponsePeriodic, ...
                thePcurrentBackgroundResponseTransient] = computeSingleConePhotocurrentResponse(...
                    eccentricityDegs, theSingleConeExcitations, temporalSupportSeconds, ...
                    nWarmUpPeriods, coneMosaicIntegrationTime, osTimeStep, pCurrentTemporalResolutionSeconds, true);
    
            % Retrieve the photocurrent response
            theConePhotocurrents(:, iConeIndex) = theConePhotoCurrentDifferentialResponse;
            theConeBackgroundPhotocurrents(iConeIndex) = theConeBackgroundPhotoCurrent;
    
            if (ismember(iConeIndex, visualizedConeIndices))
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
                disp('Hit enter to continue ...')
                pause
            end
        end % for iConeIndex
    else  % if (!debugInputConeMosaicPcurrentResponse)
    
        fprintf('Computing photocurrents\n');
        tic
        parfor iConeIndex = 1:size(theConeMosaicExcitationResponseSequence,2)
    
            if (mod(iConeIndex,100) == 0)
                % Give some feedback
                fprintf('Computing photocurrent for cone %d of %d\n', iConeIndex, size(theConeMosaicExcitationResponseSequence,2));
            end
    
            % Retrieve the cone excitation response
            theSingleConeExcitations = theConeMosaicExcitationResponseSequence(:,iConeIndex);
    
            % Compute photocurrent for this cone making the cone excitation response periodic by concatenating nWarmUpPeriods
            [~, theConePhotoCurrentDifferentialResponse, ...
                theConeBackgroundPhotoCurrent] = computeSingleConePhotocurrentResponse(...
                    eccentricityDegs, theSingleConeExcitations, temporalSupportSeconds, nWarmUpPeriods, ...
                    coneMosaicIntegrationTime, osTimeStep, pCurrentTemporalResolutionSeconds, true);
    
            % Retrieve the photocurrent response
            theConePhotocurrents(:, iConeIndex) = theConePhotoCurrentDifferentialResponse;
            theConeBackgroundPhotocurrents(iConeIndex) = theConeBackgroundPhotoCurrent;
        end % parfor iConeIndex
    
        fprintf('Photocurrents computed in %2.1f minutes\n', toc/60);
    end  % if (!debugInputConeMosaicPcurrentResponse)

end


function [temporalSupportPhotocurrent, theConePhotoCurrentDifferentialResponse, backgroundPhotocurrent, ...
    theConeExcitationsPeriodic, photocurrentResponseTimeAxis, ...
    thePcurrentDifferentialPeriodicResponse, thePcurrentBackgroundResponseTransient] = computeSingleConePhotocurrentResponse(...
    eccentricityDegs, theConeExcitations, temporalSupportSeconds, ...
    nWarmUpPeriods, coneMosaicIntegrationTime, osTimeStep, pCurrentTemporalResolutionSeconds, skipAssertions)

% Make the cone excitation response periodic by concatenating nWarmUpPeriods
theConeExcitations = theConeExcitations(1:end-1);
theConeExcitationsPeriodic = reshape(theConeExcitations, [1 numel(theConeExcitations)]);
for k = 1:nWarmUpPeriods
    theConeExcitationsPeriodic = cat(2, theConeExcitationsPeriodic, theConeExcitationsPeriodic(1:numel(theConeExcitations)));
end

% Convert excitation counts to excitation rates
theConeExcitationsRatePeriodic = theConeExcitationsPeriodic(:) / coneMosaicIntegrationTime;
backgroundConeExcitationRate = mean(theConeExcitationsRatePeriodic);

% Compute the pCurrent response to the period stimulus
[thePcurrentDifferentialPeriodicResponse, photocurrentResponseTimeAxis, thePcurrentBackgroundResponseTransient] = ...
    cMosaic.photocurrentFromConeExcitationRateUsingBiophysicalOSmodel(...
    sqrt(sum(eccentricityDegs(:).^2)), ...
    theConeExcitationsRatePeriodic, ...
    backgroundConeExcitationRate, ...
    coneMosaicIntegrationTime, ...
    pCurrentTemporalResolutionSeconds, ...
    'osTimeStepSeconds', osTimeStep, ...
    'skipAssertions', skipAssertions);

backgroundPhotocurrent = thePcurrentBackgroundResponseTransient(end);

% Only keep pCurrent response during the last stimulus period
dToriginal = temporalSupportSeconds(2)-temporalSupportSeconds(1);
dT = photocurrentResponseTimeAxis(2)-photocurrentResponseTimeAxis(1);
tOneStimulusCycle = temporalSupportSeconds(end)-temporalSupportSeconds(1);
idx = find(photocurrentResponseTimeAxis >= photocurrentResponseTimeAxis(end)-(tOneStimulusCycle+0.5*dToriginal-dT));

theConePhotoCurrentDifferentialResponse = thePcurrentDifferentialPeriodicResponse(idx);
temporalSupportPhotocurrent = photocurrentResponseTimeAxis(idx);

% Time support starts at 0 msec
temporalSupportPhotocurrent = temporalSupportPhotocurrent - temporalSupportPhotocurrent(1);
end