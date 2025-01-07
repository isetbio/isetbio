%% Initialize
ieInit;

%%
% Choose scenarios to run
%{ 
 Possible rndCodePaths for scenario: 'ConeMosaicInitialization'
    'precomputed lattice passing randomSeed';
    'custom coneData without randomSeed';
    'custom coneData with randomSeed';
    'regenerate lattice without randomSeed';
    'regenerate lattice with randomSeed';
%}

%{ 
Possible rndCodePaths for scenario: FixationalEyeMovements
    'emGenSequence without randomSeed';
    'emGenSequence passing randomSeed';
%}

%{ 
Possible rndCodePaths for scenario: ConeMosaicCompute
    'singleOInoFixationalEMrandomNoise';
    'singleOInoFixationalEMfrozenNoise';
    'singleOIwithFixationalEMrandomNoise';
%}


% Scenarios to run
scenariosList = {};

% Add scenarios to generate a list of scenarios to run
% (A) Cone mosaic initialization scenarios
%scenariosList{size(scenariosList,1)+1,1} = 'ConeMosaicInitialization';
%scenariosList{size(scenariosList,1),2} = 'custom coneData without randomSeed';

scenariosList{size(scenariosList,1)+1,1} = 'ConeMosaicInitialization';
scenariosList{size(scenariosList,1),2} = 'precomputed lattice passing randomSeed';


% (B) Fixational eye movement scenarios
%scenariosList{size(scenariosList,1)+1,1} = 'FixationalEyeMovements';
%scenariosList{size(scenariosList,1),2} = 'emGenSequence passing randomSeed';

% (C) Compute scenarios
%scenariosList{size(scenariosList,1)+1,1} = 'ConeMosaicCompute';
%scenariosList{size(scenariosList,1),2} = 'singleOInoFixationalEMrandomNoise';

%scenariosList{size(scenariosList,1)+1,1} = 'ConeMosaicCompute';
%scenariosList{size(scenariosList,1),2} = 'singleOInoFixationalEMfrozenNoise';

scenariosList{size(scenariosList,1)+1,1} = 'ConeMosaicCompute';
scenariosList{size(scenariosList,1),2} = 'singleOIwithFixationalEMrandomNoise';



% Change directory to the intercepted functions so we can run these instead
% of the original ones
theRootDir = fullfile(strrep(isetRootPath, 'isetcam', ''), 'isetbio', 'interceptedFunctions');
cd(theRootDir);


% Struct with info that we pass to intercepted rng
global rngTrackingInfo

% Generate UIfigure and save it in global variable rngTrackingInfo
hFig = uifigure();
set(hFig, 'Position', [10 10 1520 500]);
rngTrackingInfo.callingStackUIFigure = hFig;

% Run all scenarios from the scenariosList
for iScenario = 1:size(scenariosList,1)

    % Update rngTrackingInfo
    rngTrackingInfo.scenarioBeingRun = scenariosList{iScenario,1};
    rngTrackingInfo.rngCodePathToRun = scenariosList{iScenario,2};

    % Initialize the rng call number
    rngTrackingInfo.callNo = 0;

    % Run the current scenario
    switch (rngTrackingInfo.scenarioBeingRun)
        case  'ConeMosaicInitialization'
            theConeMosaic = runConeMosaicInitializationScenario(rngTrackingInfo.rngCodePathToRun);
    
        case 'FixationalEyeMovements'
            runFixationalEMgenerationScenario(rngTrackingInfo.rngCodePathToRun, theConeMosaic);
    
        case 'ConeMosaicCompute'
            runConeMosaicComputeScenario(rngTrackingInfo.rngCodePathToRun, theConeMosaic);

        otherwise
            error('Unknown scenario: ''%s''.', rngTrackingInfo.scenarioBeingRun);

    end %switch (scenario)
end % for iScenario


%
%
% SCENARIO FUNCTIONS
%
%

function runConeMosaicComputeScenario(rngCodePathToRun, theConeMosaic)

    switch (rngCodePathToRun)
        case 'singleOInoFixationalEMrandomNoise'

            % Set the noise flag to 'random'
            theConeMosaic.noiseFlag = 'random';

            % Generate test scene and OI appropriate for theConeMosaic
            [theScene, theOI] = generateConeMosaicComputeComponents(theConeMosaic, false);


            % Compute the cone mosaic response
            [theNoiseFreeAbsorptions, theNoisyAbsorptionInstances, ~, ~, temporalSupportSeconds] = ...
                theConeMosaic.compute(theOI);
            
            % Visualize everything
            visualizeConeMosaicComputeComponents(theScene, theOI, theConeMosaic, ...
                theNoiseFreeAbsorptions, theNoisyAbsorptionInstances, temporalSupportSeconds);

            % rng called from:
            %   cMosaic.compute (line 401) which calls noisyInstances with passed randomSeed
            %   cMosaic.noisyInstances (line 20) which shuffles the seed (rng('shuffle'))

            % rng called 2nd time from
            %    cMosaic.noisyInstances (line 33) which calls iePoisson
            %    iePoisson (line 109) which gets the current seed (seed = rng);

        case 'singleOInoFixationalEMfrozenNoise'

            % Set the noise flag to 'frozen'
            theConeMosaic.noiseFlag = 'frozen';

            % Generate test scene and OI appropriate for theConeMosaic
            [theScene, theOI] = generateConeMosaicComputeComponents(theConeMosaic, false);

            % Compute the cone mosaic response
            [theNoiseFreeAbsorptions, theNoisyAbsorptionInstances, ~, ~, temporalSupportSeconds] = ...
                theConeMosaic.compute(theOI);
            
            % Visualize everything
            visualizeConeMosaicComputeComponents(theScene, theOI, theConeMosaic, ...
                theNoiseFreeAbsorptions, theNoisyAbsorptionInstances, temporalSupportSeconds);

            % rng called from:
            %   cMosaic.compute (line 401) which calls noisyInstances with passed randomSeed
            %   cMosaic.noisyInstances (line 18) which sets the passed seed (rng(seed));

            % rng called 2nd time from
            %    cMosaic.noisyInstances (line 33) which calls iePoisson
            %    iePoisson (line 107) which sets the passed seed (rng(seed));


        case 'singleOIwithFixationalEMrandomNoise'

            % Set the noise flag to 'random'
            theConeMosaic.noiseFlag = 'random';

            % Generate test scene and OI appropriate for theConeMosaic
            [theScene, theOI, fixationalEMobj] = generateConeMosaicComputeComponents(theConeMosaic, true);


            % Compute the cone mosaic response
            [theNoiseFreeAbsorptions, theNoisyAbsorptionInstances, ~, ~, temporalSupportSeconds] = ...
                theConeMosaic.compute(theOI, ...
                'withFixationalEyeMovements', true);
            
            % Visualize everything
            visualizeConeMosaicComputeComponents(theScene, theOI, theConeMosaic, ...
                theNoiseFreeAbsorptions, theNoisyAbsorptionInstances, temporalSupportSeconds);

            % rng called from:
            %   fixationalEMObj.computeForCmosaic (line 66) which calls
            %   fixationalEMObj.compute (line 65) which shuffles the seed (rng('shuffle'))

            % rng called 2nd time from:
            %   cMosaic.compute (line 401) which calls noisyInstances with passed randomSeed
            %   cMosaic.noisyInstances (line 20) which shuffles the seed (rng('shuffle'))

            % rng called 3rd time from
            %    cMosaic.noisyInstances (line 33) which calls iePoisson
            %    iePoisson (line 109) which gets the current seed (seed = rng);


        otherwise
            error('Unknown case to run: ''%s''.', rngCodePathToRun)
    end  % switch (rngCodePathToRun)


end

function runFixationalEMgenerationScenario(rngCodePathToRun, theConeMosaic)

    eyeMovementDurationSeconds = 0.2;
    switch (rngCodePathToRun)
        case 'emGenSequence without randomSeed'
            theConeMosaic.emGenSequence(eyeMovementDurationSeconds, ...
                'microsaccadeType', 'none', ...
                'centerPaths', true, ...
                'nTrials', 2);

            % rng called from:
            %   cMosaic.emGenSequence (line 52) to save the current rng state struct

            % rng called 2nd time from:
            %   cMosaic.emGenSequence (line 105) which calls obj.fixEMobj.computeForCmosaic().compute
            %   obj.fixEMobj.computeForCmosaic.compute (line 66) sets the passed random seed
            
            % rng called 3rd time from:
            %   cMosaic.emGenSequence (line 114) to restore the previous rng state struct
            
            
        case 'emGenSequence passing randomSeed'
            theConeMosaic.emGenSequence(eyeMovementDurationSeconds, ...
                'microsaccadeType', 'none', ...
                'centerPaths', true, ...
                'nTrials', 2, ...
                'randomSeed', 2345);

            % rng called from:
            %   cMosaic.emGenSequence (line 52) to save the current rng state struct

            % rng called 2nd time from:
            %   cMosaic.emGenSequence (line 105) which calls obj.fixEMobj.computeForCmosaic().compute
            %   obj.fixEMobj.computeForCmosaic.compute (line 68) shuffles the seed (rng('shuffle')
            
            % rng called 3rd time from:
            %   cMosaic.emGenSequence (line 114) to restore the previous rng state struct



          otherwise
            error('Unknown case to run: ''%s''.', rngCodePathToRun)
    end  % switch (rngCodePathToRun)

end




function theConeMosaic = runConeMosaicInitializationScenario(rngCodePathToRun)

    coneMosaicIntegrationTimeMilliSeconds = 1;

    switch (rngCodePathToRun)
        case 'precomputed lattice passing randomSeed'
            % CASE 1: Generate a @cMosaic passing ('randomSeed', val) key-value pair
            theConeMosaic = cMosaic(...
                'sizeDegs', [0.5 0.5], ...      
                'eccentricityDegs', [0 0], ... 
                'integrationTime', coneMosaicIntegrationTimeMilliSeconds/1000, ...    
                'randomSeed', 22 ...
                );
            % rng called from:
            %   cMosaic.cMosaic (line 556)
    
        case 'custom coneData without randomSeed'
            % CASE 2: Generate a @cMosaic from completely custom cone data WITHOUT
            % specifying a randomSeed
            [theCustomConeDataStruct, theCustomRetinalMagnification] = helperGenerateCustomConeData();
            theConeMosaic = cMosaic(...
                     'coneData', theCustomConeDataStruct, ...
                     'micronsPerDegree', theCustomRetinalMagnification, ...
                     'integrationTime', coneMosaicIntegrationTimeMilliSeconds/1000);
            
            % rng called from:
            %   cMosaic.cMosaic (line 576)
    
        case 'custom coneData with randomSeed'
            % CASE 3: Generate a @cMosaic from completely custom cone data AND
            % a specified randomSeed
            [theCustomConeDataStruct, theCustomRetinalMagnification] = helperGenerateCustomConeData();
            theConeMosaic = cMosaic(...
                     'coneData', theCustomConeDataStruct, ...
                     'micronsPerDegree', theCustomRetinalMagnification, ...
                     'integrationTime', coneMosaicIntegrationTimeMilliSeconds/1000, ... 
                     'randomSeed', 123);
            
            % rng called from:
            %   cMosaic.cMosaic (line 582)
    
        
        case 'regenerate lattice without randomSeed'
            % CASE 4: Generate a @cMosaic by regenerating the lattice (not importing it) WITHOUT
            % specifying a randomSeed
            theConeMosaic = cMosaic(...
                     'coneData', [], ...
                     'computemeshfromscratch', true, ...
                     'integrationTime', coneMosaicIntegrationTimeMilliSeconds/1000);
    
            % rng called with seed set to 1 from:
            %   cMosaic.regenerationConePositions (line 540) -> 
            %   retinalattice.generatePatch ->  (calls retinalattice.configure, which sets the randomSeed to 1)
            %   retinalattice.initialize.downSampleInitialRFpositions (line17)
    
        case 'regenerate lattice with randomSeed'
            % CASE 5: Generate a @cMosaic by regenerating the lattice (not importing it) AND
            % a specified randomSeed
            theConeMosaic = cMosaic(...
                     'coneData', [], ...
                     'computemeshfromscratch', true, ...
                     'randomSeed', 555, ...
                     'integrationTime', coneMosaicIntegrationTimeMilliSeconds/1000);
    
            % rng called with seed set to 1 from:
            %   cMosaic.regenerationConePositions (line 540) ->
            %   retinalattice.generatePatch ->  (calls retinalattice.configure, which sets the randomSeed to 1, but then gets overriden to 555)
            %   retinalattice.initialize.downSampleInitialRFpositions (line17)
    
            % rng called 2nd time with seed 555 from:
            %   cMosaic.cMosaic (line 556)
    
        otherwise
            error('Unknown case to run: ''%s''.', rngCodePathToRun)
    end % switch (rngCodePathToRun)

end % function runConeMosaicInitializationScenario()




%
% --- HELPER FUNCTIONS ----
%

function [theCustomConeDataStruct, theCustomRetinalMagnification] = helperGenerateCustomConeData()
    % All specifications are in units of retinal microns
    thePositionUnits = 'microns';

    % Specify a custom retinal magnification, here 250 microns/deg
    theCustomRetinalMagnification = 250;

    % Specify the (x,y) positions of all cones (in microns)
    % Here we are laying cones along a spiral path with cone aperture increasing
    % with distance from the center.
    angles = 0:15:3000;
    minConeDiameter = 2;
    for iCone = 1:numel(angles)
        iAngle = angles(iCone);
        radius = minConeDiameter + iCone*0.7;
        thePositions(iCone,:) = radius * [ cosd(iAngle) sind(iAngle)];
        theConeApertureDiameters(iCone) = minConeDiameter + (iCone*minConeDiameter*0.01)^1.2;
        switch (mod(iCone,11))
            case {1,2,4,5,7,8,9}
                theTypes(iCone) = cMosaic.LCONE_ID;  % an L-cone
            case {0,3,6}
                theTypes(iCone) = cMosaic.MCONE_ID;  % an M-cone
            case 10
                theTypes(iCone) = cMosaic.SCONE_ID;  % an S-cone
        end
    end
    

    % Generate struct with custom cone data
    theCustomConeDataStruct = struct(...
        'positionUnits', thePositionUnits, ...
        'positions', thePositions, ...
        'types', theTypes,...
        'lightGatheringApertureDiameters', theConeApertureDiameters ...                      
        );

end


function [theScene, theOI, fixationalEMObj] = generateConeMosaicComputeComponents(theConeMosaic, generateFixationalEMobj)

    vParams = vernierP;
    vParams.barWidth = 2;
    vParams.offset = 3;
    theScene = sceneCreate('vernier', vParams.sceneSz(1), vParams.barWidth, vParams.offset);
    theScene = sceneSet(theScene, 'fov', 0.3);

    % Best Polans subject
    opticsDataBase = 'Polans2015';
    subjectRankOrder = 1;
    rankedSujectIDs = PolansOptics.constants.subjectRanking;
    testSubjectID = rankedSujectIDs(subjectRankOrder);

    % Determine if we need to subtract the subject's central refraction
    subtractCentralRefraction = PolansOptics.constants.subjectRequiresCentralRefractionCorrection(testSubjectID);

    % Generate theOI
    oiEnsemble = theConeMosaic.oiEnsembleGenerate(theConeMosaic.eccentricityDegs, ...
        'zernikeDataBase', opticsDataBase, ...
        'subjectID', testSubjectID, ...
        'pupilDiameterMM', 3.0, ...
        'zeroCenterPSF', true, ...
        'subtractCentralRefraction', subtractCentralRefraction);
    theOI = oiEnsemble{1,1};

    % Compute the retinal image
    theOI = oiCompute(theOI, theScene);

    if (generateFixationalEMobj)
        % Initialize
        fixationalEMObj = fixationalEM;              % Instantiate a fixationalEM object
        fixationalEMObj.microSaccadeType = 'none';   % No microsaccades, just drift
    
        % Compute number of eye movements
        nTrials = 1;
        frameDurationSeconds = theConeMosaic.integrationTime;
        stimDurationSeconds = frameDurationSeconds * 30;
        eyeMovementsPerTrial = stimDurationSeconds/frameDurationSeconds;

        % Generate the em sequence for the passed cone mosaic,
        % which results in a time step equal to the integration time of theConeMosaic
        fixationalEMObj.computeForCmosaic(...
            theConeMosaic, eyeMovementsPerTrial,...
            'nTrials' , nTrials);

        % Set the fixational eye movements into the cone mosaic
        theConeMosaic.emSetFixationalEMObj(fixationalEMObj);

    else
        fixationalEMObj = [];
    end

end


function visualizeConeMosaicComputeComponents(theScene, theOI, theConeMosaic, ...
    theNoiseFreeAbsorptions, theNoisyAbsorptionInstances, temporalSupportSeconds)

    activationRange = [min(theNoiseFreeAbsorptions(:)) max(theNoiseFreeAbsorptions(:))];

    scenePixelWidthDegs = sceneGet(theScene, 'wangular resolution');
    sceneWidthPixels = sceneGet(theScene, 'cols');
    sceneSupport = (1:sceneWidthPixels)*scenePixelWidthDegs;
    sceneSupport = sceneSupport - mean(sceneSupport);

    oiPixelWidthDegs = oiGet(theOI, 'wangular resolution');
    oiWidthPixels = oiGet(theOI, 'cols');
    oiSupport = (1:oiWidthPixels)*oiPixelWidthDegs;
    oiSupport = oiSupport - mean(oiSupport);

    hFig = figure(1); clf;
    set(hFig, 'Position', [10 10 1000 1000]);

    % The stimulus
    ax1 = subplot(2,2,1);
    image(ax1, sceneSupport, sceneSupport, sceneGet(theScene, 'rgbimage'));
    axis(ax1, 'image');
    set(ax1, 'Color', [0 0 0], 'XLim', 0.5*theConeMosaic.sizeDegs(1)*[-1 1], 'YLim', 0.5*theConeMosaic.sizeDegs(2)*[-1 1] );
    set(ax1, 'FontSize', 16);
    title(ax1, 'stimulus');

    % The optical image
    ax2 = subplot(2,2,2);
    image(ax2, oiSupport, oiSupport, oiGet(theOI, 'rgbimage'));
    axis(ax2, 'image');
    set(ax2, 'Color', [0 0 0], 'XLim', 0.5*theConeMosaic.sizeDegs(1)*[-1 1], 'YLim', 0.5*theConeMosaic.sizeDegs(2)*[-1 1] );
    title(ax2, 'retinal image');
    set(ax2, 'FontSize', 16);

    
    % The cone mosaic response
    iTrial = 1; iTimePoint = 1;
    [nTrials, nTimePoints, nCones] = size(theNoiseFreeAbsorptions);

    if (nTimePoints > 1)
        for iTimePoint = 1:numel( temporalSupportSeconds)
            ax3 = subplot(2,2,3);
            theConeMosaic.visualize(...
                'figureHandle', hFig, ...
                'axesHandle', ax3, ...
                'activation', theNoiseFreeAbsorptions(iTrial, iTimePoint,:), ...
                'activationRange', activationRange, ...
                'plotTitle', sprintf('noise-free response (%2.1f msec)', temporalSupportSeconds(iTimePoint)*1000), ...
                'fontSize', 16 ...
                );

            ax4 = subplot(2,2,4);
            theConeMosaic.visualize(...
                'figureHandle', hFig, ...
                'axesHandle', ax4, ...
                'activation', theNoisyAbsorptionInstances(iTrial, iTimePoint,:), ...
                'activationRange', activationRange, ...
                'plotTitle', sprintf('noisy response instance (%2.1f msec)', temporalSupportSeconds(iTimePoint)*1000), ...
                'fontSize', 16 ...
                );

            drawnow;
        end

    else
        ax3 = subplot(2,2,3);
        theConeMosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax3, ...
            'activation', theNoiseFreeAbsorptions(iTrial, iTimePoint,:), ...
            'activationRange', activationRange, ...
            'plotTitle', 'noise-free response', ...
            'fontSize', 16 ...
            );
        ax4 = subplot(2,2,4);
        theConeMosaic.visualize(...
            'figureHandle', hFig, ...
            'axesHandle', ax4, ...
            'activation', theNoisyAbsorptionInstances(iTrial, iTimePoint,:), ...
            'activationRange', activationRange, ...
            'plotTitle', 'noisy response instance', ...
            'fontSize', 16 ...
            );
    end


    
end
