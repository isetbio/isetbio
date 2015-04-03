function varargout = v_DisplayColorConversion(varargin)
%
% Validate display calibration color conversion against PTB.
%
% In the graphs, the comparison is ballpark good and we are happy enough with that.
% We've saved out what the two do, after comparing the graphs here.
%
% See also v_IrradianceIsomerizations

varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

%% Initialize ISET
s_initISET;

%% Some informative text
UnitTest.validationRecord('SIMPLE_MESSAGE', 'Compare isetbio and PTB color conversion.');

%% Remove the Brainard lab PTB overrides folder from the path
%
% This prevents this code from using the new BL object oriented PTB
% overrides.  We could include these inside of isetbio, but the risk
% is that whether this program worked or not would depend on whether
% isetbio was before or after PTB on the user's path, something we
% don't want to have to deal with.  Elsewhere (not in isetbio), we have
% established that the PTB and BrainardLab routines do the same thing.
[removedFolderFromCurrentPath, originalPath] = removeBrainardLabPTBOverrides();

try
    %% Overview
    %
    % What we are going to do in this script is to create
    % uniform field corresponding to specified display and
    % rgb values and pump it through isetbio to get cone 
    % isomerization rates.  Then we are going to do the
    % calculation for the same display, cone fundamentals, and
    % rgb values in PTB.  We start with the isetbio calculation.
    
    %% Create isetbio display 
    displayToTest = 'LCD-Apple';
    d = displayCreate(displayToTest);
    temp = displayGet(d, 'gamma table');
    if (size(temp,2) ~=  3)
        error('Don''t feel like dealing with a display with other than 3 primaries for this test.');
    end
    
    %% Create a uniform field radiance image in ISETBIO, using the display and the
    % passed RGB values.
    RGBToTest = [0.2 0.73 0.42]';
    theRGBImage = ones(20,20,3);
    for i1 = 1:3
        theRGBImage(:,:,i1) = round(255*RGBToTest(i1));
    end
    sceneDegrees = 20;  % need large field
    scene = sceneFromFile(theRGBImage,'rgb',[],d);
    scene = sceneSet(scene,'fov', sceneDegrees);
    if (runTimeParams.generatePlots)
        vcAddAndSelectObject(scene); sceneWindow;
    end
    
    %% Compute optial image, for delta function optics.
    %
    % There may be a better way to specify delta function
    % optics than what I do here, but the kluge I'm using
    % here of setting the OTF to all ones this works.
    oi = oiCreate('human');
    optics = oiGet(oi,'optics');
    OTFData = opticsGet(optics,'otfdata');
    OTFDeltaData = ones(size(OTFData));
    optics = opticsSet(optics,'otfdata',OTFDeltaData);
    oi = oiSet(oi,'optics',optics);
    oi = oiCompute(scene,oi);
    if (runTimeParams.generatePlots)
        vcAddAndSelectObject(oi); oiWindow;
    end
    
    %% Compute sensor image
    sensorDegrees = 4;
    sensorROIPixels = 10;
    sensor = sensorCreate('human');
    sensor = sensorSet(sensor, 'noise flag', 0);
    sensor = sensorSet(sensor,'exp time',1);
    sensor = sensorSet(sensor,'rows',128);
    sensor = sensorSet(sensor,'cols',128);
    [sensor, ~] = sensorSetSizeToFOV(sensor,sensorDegrees,scene,oi);
    sensor = sensorCompute(sensor,oi);
    if (runTimeParams.generatePlots)
        vcAddAndSelectObject(sensor); sensorWindow('scale',1);
    end
    
    %% Get the LMS isomerization rates out of the sensor image
    %
    % Pull out responses of each cone type within ROI. I am doing this by
    % brute force, because I can't find quite the right combination of ROI
    % gets from the sensor image.
    isetbioIsomerizationsArray = sensorGet(sensor,'photons');    
    sensorCFA = sensorGet(sensor,'cfa');
    sensorSizePixels = sensorGet(sensor,'size');
    rect = round([sensorSizePixels(2)/2,sensorSizePixels(1)/2,sensorROIPixels,sensorROIPixels]);
    sensorRoiLocs = ieRoi2Locs(rect);
    nLocs = size(sensorRoiLocs,1);
    sumIsomerizations = zeros(3,1);
    nSummed = zeros(3,1);    
    for jj = 1:nLocs
        % A type of 1 in the CFA is blank, so we subtract 1 from the number
        % in the CFA and skip any that end up as 0
        coneType = sensorCFA.pattern(sensorRoiLocs(jj,1),sensorRoiLocs(jj,2))-1;
        if (coneType > 0)
            sumIsomerizations(coneType) = sumIsomerizations(coneType)+isetbioIsomerizationsArray(sensorRoiLocs(jj,1),sensorRoiLocs(jj,2));
            nSummed(coneType) = nSummed(coneType) + 1;
        end
    end
    isetbioLMSIsomerizations = sumIsomerizations ./ nSummed;
    
    %% Now we're going to get the isomerizations the PTB way
    % Get cone fundamentals out of isetbio sensor, and convert
    % to energy units.  We'll use these for PTB calcs
    S_cones = WlsToS(sensorGet(sensor,'wave'));
    T_conesQE = sensorGet(sensor,'spectral qe')';
    T_conesQE = T_conesQE(2:4,:);
    T_cones = EnergyToQuanta(S_cones,T_conesQE')';

    %% Create PTB calibration structure
    gammaTable = displayGet(d, 'gamma table');
    nInputLevels = size(gammaTable,1);
    wave = displayGet(d, 'wave');
    spd  = displayGet(d, 'spd primaries');
    originalGammaTableLength = size(gammaTable,1);
    originalSettingsValues   = linspace(0,1,originalGammaTableLength);
    dotsPerMeter = displayGet(d, 'dots per meter');
    screenSizeInPixels = [1920 1080];
    PTBcal = ptb.GeneratePsychToolboxCalStruct(...
        'name', displayGet(d, 'name'), ...
        'gammaTable', gammaTable, ...
        'wave', wave, ...
        'spd', spd, ...
        'screenSizeInPixels', screenSizeInPixels, ...
        'dotsPerMeter', dotsPerMeter ...
        );
    
    % Kluge.  Should set this properly in creation routine
    PTBcal.P_ambient = zeros(length(wave),1);
    
    %% Initialize PTB cal
    PTBcal = CalibrateFitGamma(PTBcal, nInputLevels);
    gammaMethod = 1;
    PTBcal = SetGammaMethod(PTBcal, gammaMethod, nInputLevels);
    PTBcal = SetSensorColorSpace(PTBcal,T_cones,S_cones);
    
    %% PTB conversion to isomerization rates
    ptbLMSIsomerizations = SettingsToSensor(PTBcal,RGBToTest);
    
    %% NOW NEED TO WORK THROUGH AND GET UNITS RIGHT IN THE TWO CALCULATIONS
    % AND FIGURE OUT WHY RATIO IS NOT QUITE CONSTANT.
    ptbLMSIsomerizations./isetbioLMSIsomerizations
    
    % Log the data
    dataStruct = struct( ...
        'displayName', displaysToTest{displayIndex}, ...
        'originalSettingsValues', originalSettingsValues, ...
        'gammaTable', gammaTable, ...
        'nInputLevels', nInputLevels, ...
        'settingsValues', settingsValues, ...
        'PTBinverseGamma', PTBinverseGamma, ...
        'ISETBIOinverseGamma', ISETBIOinverseGamma ...
        );
    
catch err
    % Restore original path and rethrow error
    if (~isempty(removedFolderFromCurrentPath))
        path(originalPath);
    end
    rethrow(err);
end

%% Restore original path
if (~isempty(removedFolderFromCurrentPath))
    path(originalPath);
end

%% Save validation operations
UnitTest.validationData('inversionData', dataStruct);

%% Plot
if (runTimeParams.generatePlots)
    h = figure(1);
    set(h, 'Position', [10 10 1704 1196]);
    clf;
    
    subplotIndex = 0;
    subplotIndex = subplotIndex + 1;
    subplot(numel(displaysToTest),numel(gammaTableLengthsToTest),subplotIndex);
    hold on;
    
    % Plot results for RED channel only
    channelIndex = 1;
    
    % original gamma
    plot(dataStruct.originalSettingsValues, dataStruct.gammaTable(:,channelIndex), 'g.');
    
    % PTB-inverted
    plot(dataStruct.settingsValues, dataStruct.PTBinverseGamma(:,channelIndex), 'r.');
    
    % ISETBIO-inverted
    plot(dataStruct.settingsValues, dataStruct.ISETBIOinverseGamma(:,channelIndex), 'b.');
    
    % ISETBIO-inverted
    plot(dataStruct.settingsValues, 0.5+dataStruct.PTBinverseGamma(:,channelIndex)-dataStruct.ISETBIOinverseGamma(:,channelIndex), 'k-', 'LineWidth', 2.0);
    
    hold off;
    
    h_legend = legend('Orig. LUT', 'PTB inverse', 'ISETBIO inverse', '0.5+(PTB-ISETBIO)', 'Location', 'NorthWest');
    set(h_legend,'FontName', 'Menlo', 'FontSize',12);
    
    set(gca, 'XLim', [0 1], 'YLim', [0 1], 'XTick', [0:0.1:1.0], 'YTick', [0:0.1:1.0], 'FontSize', 12, 'FontName', 'Helvetica');
    xlabel('gamma in','FontSize', 12, 'FontName', 'Helvetica', 'FontWeight', 'bold');
    ylabel('gamma out', 'FontSize', 12, 'FontName', 'Helvetica', 'FontWeight', 'bold');
    axis 'square'; box on; grid on;
    title(sprintf('%s - (%d levels)', dataStruct.displayName, dataStruct.nInputLevels), 'FontSize', 16, 'FontName', 'Helvetica', 'FontWeight', 'bold');
    hold off
end
end

% Helper method to remove the BrainardLabPTBOverrides folder if it exists on the current path.
% We want to remove this override, so we can use the original PTB functions
% without the CalStructOBJ (which is found only on our computers).
function [removedFolderFromCurrentPath, originalPath] = removeBrainardLabPTBOverrides

removedFolderFromCurrentPath = '';
originalPath = path;

% Folder to remove the Overrides/PTB-3 folder from the current path, if the folder exists on the path
PTBoverridesDirToRemoveFromPath = '/Users/Shared/Matlab/Toolboxes/BrainardLabToolbox/Overrides/PTB-3';

% determine if the PTBoverridesDirToRemoveFromPath is in the current path
pathCell = regexp(path, pathsep, 'split');
onPath = any(strcmpi(PTBoverridesDirToRemoveFromPath, pathCell));
if (onPath)
    rmpath(PTBoverridesDirToRemoveFromPath);
    removedFolderFromCurrentPath = PTBoverridesDirToRemoveFromPath;
end
end

