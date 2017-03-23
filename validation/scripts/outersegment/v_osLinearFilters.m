function varargout = v_osLinearFilters(varargin)
%
% Compute linearFilters at different os.timeSteps and luminance levels
%
% NPC, ISETBIO TEAM, 2016

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)
%% Init
ieInit;

% Define the time axis for the simulation
stimulusSamplingInterval = 10/1000;            % 5 milliseconds
oiTimeAxis = 0:stimulusSamplingInterval:0.1;   % 0.1 seconds
oiTimeAxis = oiTimeAxis - mean(oiTimeAxis);

FOV = 1.0;
% Valid range of background luminances - result in 500 - 20,0000 isomerizations/cone/second
% regime where photocurrent model was developed on

backgroundLuminances = [3.3]; 
osTimeSteps = [0.1]/1000;

% Compute the stimulus modulation function
stimulusRampTau = 180/1000;
modulationGain = 1.0;
modulationFunction = modulationGain * exp(-0.5*(oiTimeAxis/stimulusRampTau).^2);

% Generate optics
noOptics = false;
theOI = oiGenerate(noOptics);

% Generate a cone mosaic with 1 L-, 1 M-, and 1 S-cone only
mosaicSize =  nan;                    
integrationTime = 2.5/1000;            
photonNoise = 'none';                                        
osNoise = 'none';                        
theConeMosaic = coneMosaicGenerate(mosaicSize, photonNoise, osNoise, integrationTime, osTimeSteps(1));

for iLum = 1:numel(backgroundLuminances)
    fprintf('Computing os linear filters for background luminance %2.1f cd/m2  [%d/%d]\n',  backgroundLuminances(iLum), iLum, numel(backgroundLuminances));
    theScene = uniformFieldSceneCreate(FOV, backgroundLuminances(iLum));
    % Compute the optical images
    oiBackground = oiCompute(theOI, theScene);
    oiModulated  = oiBackground;

    % Generate the sequence of optical images
    theOIsequence = oiSequence(oiBackground, oiModulated, oiTimeAxis, modulationFunction);

    % Generate eye movement path
    eyeMovementsNum = theOIsequence.maxEyeMovementsNumGivenIntegrationTime(theConeMosaic.integrationTime);
    instancesNum = 1;
    emPaths = zeros(instancesNum, eyeMovementsNum,2);

    for iStepIndex = 1:numel(osTimeSteps)
        theConeMosaic.os.timeStep = osTimeSteps(iStepIndex);
        % Compute all instances 
        [abosptionCountsPerIntegrationTime, ~, LMSfilters{iLum, iStepIndex}] = ...
            theConeMosaic.computeForOISequence(theOIsequence, ...
            'emPaths', emPaths, ...
            'currentFlag', true);

         isomerizationsPerSecond = max(abosptionCountsPerIntegrationTime(:)) / theConeMosaic.integrationTime;
         [backgroundLuminances(iLum) isomerizationsPerSecond 500 20000]
    end % iStepIndex
end

filterTimeAxis = theConeMosaic.interpFilterTimeAxis;

% Unit test validation data
UnitTest.validationRecord('SIMPLE_MESSAGE', '***** v_osLinearFilters *****');
UnitTest.validationData('LMSfilters', LMSfilters, ...
    'UsingTheFollowingVariableTolerancePairs', ...
    'LMSFilters', 1e-4);
UnitTest.validationData('filterTimeAxis', filterTimeAxis);
UnitTest.validationData('backgroundLuminances', backgroundLuminances);
UnitTest.validationData('osTimeSteps', osTimeSteps);

    
if (runTimeParams.generatePlots)
    plotIRs(osTimeSteps, backgroundLuminances, LMSfilters, filterTimeAxis)
end
end

function plotIRs(osTimeSteps, backgroundLuminances, LMSfilters, filterTimeAxis)
       
subplotPosVectors = NicePlot.getSubPlotPosVectors(...
           'rowsNum', 2, ...
           'colsNum', 3, ...
           'heightMargin',   0.04, ...
           'widthMargin',    0.05, ...
           'leftMargin',     0.05, ...
           'rightMargin',    0.00, ...
           'bottomMargin',   0.04, ...
           'topMargin',      0.04);
       
coneNames = {'L', 'M', 'S'};

hFig = figure(1); clf;
set(hFig, 'Position', [10 10 1380 1010], 'Color', [1 1 1]);

colorIR = (jet(numel(osTimeSteps))).^2;
for iLum = 1:numel(backgroundLuminances)
    legends = {};
    for iStepIndex = 1:numel(osTimeSteps)
        legends{numel(legends)+1} = sprintf('os.timeStep: %2.3fms', osTimeSteps(iStepIndex)*1000);
        color = squeeze(colorIR(iStepIndex,:));
        for coneIndex = 1:3  
            subplot('Position', subplotPosVectors(1,coneIndex).v);
            normFilter = squeeze(LMSfilters{iLum, iStepIndex}(:,coneIndex));
            normFilter = normFilter / max(normFilter);
            plot(filterTimeAxis, normFilter, 'k-', 'Color', color, 'LineWidth', 1.5);
            if (iStepIndex == 1) 
                hold on;
            end
            if (iStepIndex == numel(osTimeSteps))
                hold off;
            end
            set(gca, 'XLim', [filterTimeAxis(1) filterTimeAxis(end)], 'YLim', [-0.1 1.05], 'YTick', (0:0.1:1.0), 'XTick', (0:0.05:0.5), 'FontSize', 12);
            grid on; box on;
            hL = legend(legends);
            set(hL, 'FontSize', 12);
            title(sprintf('%s-cone, os.timeStep: %2.3fms, lum: %2.1f cd/m2', coneNames{coneIndex}, osTimeSteps(iStepIndex)*1000, backgroundLuminances(iLum)), ...
                'FontSize', 14, 'FontWeight', 'bold');

            subplot('Position', subplotPosVectors(2,coneIndex).v);
            filter = squeeze(LMSfilters{iLum, iStepIndex}(:,coneIndex));
            plot(filterTimeAxis, filter, 'k-', 'Color', color, 'LineWidth', 1.5);
            if (iStepIndex == 1) 
                hold on;
            end
            if (iStepIndex == numel(osTimeSteps))
                hold off;
            end
            set(gca, 'XLim', [filterTimeAxis(1) filterTimeAxis(end)], 'YLim', [-0.1 1.0]*0.18, 'YTick', (0:0.1:1.0)*0.2, 'XTick', (0:0.05:0.5), 'FontSize', 12);
            grid on; box on;
            hL = legend(legends);
            set(hL, 'FontSize', 12);
            title(sprintf('%s-cone, os.timeStep: %2.3fms, lum: %2.1f cd/m2', coneNames{coneIndex}, osTimeSteps(iStepIndex)*1000, backgroundLuminances(iLum)), ...
                'FontSize', 14, 'FontWeight', 'bold');
        end % coneIndex
        drawnow
    end % iStepIndex
end %iLum
end

function theConeMosaic = coneMosaicGenerate(mosaicSize, photonNoise, osNoise, integrationTime, osTimeStep)
    % Default human mosaic
    theConeMosaic = coneMosaic;
    
    % Adjust size
    if isnan(mosaicSize)
        % Generate a human cone mosaic with 1L, 1M and 1S cone
        theConeMosaic.rows = 1;
        theConeMosaic.cols = 3;
        theConeMosaic.pattern = [2 3 4];
    else
        theConeMosaic.setSizeToFOV(mosaicSize);
    end
    
    % Set the noise
    theConeMosaic.noiseFlag = photonNoise;

    % Set the integrationTime
    theConeMosaic.integrationTime = integrationTime;
    
    % Generate the outer-segment object to be used by the coneMosaic
    theOuterSegment = osLinear();
    theOuterSegment.noiseFlag = osNoise;
    
    % Set a custom timeStep, for @osLinear we do not need the default 0.1 msec
    theOuterSegment.timeStep = osTimeStep;

    % Couple the outersegment object to the cone mosaic object
    theConeMosaic.os = theOuterSegment;
end


function theOI = oiGenerate(noOptics)
    % Generate optics
    if (noOptics)
        theOI = oiCreate('diffraction limited');
        optics = oiGet(theOI,'optics');           
        optics = opticsSet(optics,'fnumber',0);
        optics = opticsSet(optics, 'off axis method', 'skip');
        theOI = oiSet(theOI,'optics', optics);
    else
        theOI = oiCreate('human');
    end
end

function uniformScene = uniformFieldSceneCreate(FOV, meanLuminance)
    uniformScene = sceneCreate('uniform equal photon', 128);
    % square scene with desired FOV
    uniformScene = sceneSet(uniformScene, 'wAngular', FOV);
    % 1 meter away
    uniformScene = sceneSet(uniformScene, 'distance', 1.0);
    % adjust radiance according to desired  mean luminance
    uniformScene = sceneAdjustLuminance(uniformScene, meanLuminance);
end

