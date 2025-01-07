%% Initialize
ieInit;
theRootDir = fullfile(strrep(isetRootPath, 'isetcam', ''), 'isetbio');
cd(theRootDir);

    
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


% Scenarios to run
% (1) Cone mosaic initialization scenarios
scenarioList{1,1} = 'ConeMosaicInitialization';
scenarioList{1,2} = 'precomputed lattice passing randomSeed';

% (2) Fixational eye movement scenarios
scenarioList{size(scenarioList,1)+1,1} = 'FixationalEyeMovements';
scenarioList{size(scenarioList,1),2} = 'emGenSequence passing randomSeed';

% (3) Compute scenarios
%scenarioList{size(scenarioList,1)+1,1} = 'ConeMosaicCompute';
%scenarioList{size(scenarioList,1),2} = '??';


% Struct with info on what is being currently run that we pass to intercepted rng
global rngTrackingInfo

hFig = uifigure();
set(hFig, 'Position', [10 10 1520 400]);
rngTrackingInfo.callingStackUIFigure = hFig;

for iScenario = 1:size(scenarioList,1)

    % Update rngTrackingInfo
    rngTrackingInfo.scenarioBeingRun = scenarioList{iScenario,1};
    rngTrackingInfo.rngCodePathToRun = scenarioList{iScenario,2};
    rngTrackingInfo.callNo = 0;

    switch (rngTrackingInfo.scenarioBeingRun)
        case  'ConeMosaicInitialization'
            theConeMosaic = runConeMosaicInitializationScenario(rngTrackingInfo.rngCodePathToRun);
    
        case 'FixationalEyeMovements'
            runFixationalEMgenerationScenario(rngTrackingInfo.rngCodePathToRun, theConeMosaic);
    
        otherwise
            error('Unknown scenario: ''%s''.', rngTrackingInfo.scenarioBeingRun);

    end %switch (scenario)
end % for iScenario


%
%
% SCENARIO FUNCTIONS
%
%

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

    switch (rngCodePathToRun)
        case 'precomputed lattice passing randomSeed'
            % CASE 1: Generate a @cMosaic passing ('randomSeed', val) key-value pair
            theConeMosaic = cMosaic(...
                'sizeDegs', [0.5 0.5], ...      
                'eccentricityDegs', [0 0], ... 
                'integrationTime', 10/1000, ...    
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
                     'micronsPerDegree', theCustomRetinalMagnification);
            
            % rng called from:
            %   cMosaic.cMosaic (line 576)
    
        case 'custom coneData with randomSeed'
            % CASE 3: Generate a @cMosaic from completely custom cone data AND
            % a specified randomSeed
            [theCustomConeDataStruct, theCustomRetinalMagnification] = helperGenerateCustomConeData();
            theConeMosaic = cMosaic(...
                     'coneData', theCustomConeDataStruct, ...
                     'micronsPerDegree', theCustomRetinalMagnification, ...
                     'randomSeed', 123);
            
            % rng called from:
            %   cMosaic.cMosaic (line 582)
    
        
        case 'regenerate lattice without randomSeed'
            % CASE 4: Generate a @cMosaic by regenerating the lattice (not importing it) WITHOUT
            % specifying a randomSeed
            theConeMosaic = cMosaic(...
                     'coneData', [], ...
                     'computemeshfromscratch', true);
    
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
                     'randomSeed', 555);
    
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

