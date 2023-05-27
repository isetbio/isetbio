function varargout = v_PTBcalStructToIsetbioDisplayObjectAndBack(varargin)
%
% Validate conversion of PTBcalStruct to isetbio display object and back.
% This compares the PTB reconstruction after conversion into isetbio and
% back with the originally read PTB cal stucture.
%
% See also v_IrradianceIsomerizations, v_DisplayColorConversion

    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Function implementing the isetbio validation code
function ValidationFunction(runTimeParams)

    %% Initialize ISET
    close all; ieInit;
    
    %% Some informative text
    UnitTest.validationRecord('SIMPLE_MESSAGE', 'Conversion of PTB cal struct to isetbio display object.');
    
    %% Load PTB calibration data for the SONY PVM2541 display
    calDir = fullfile(isetbioDataPath, 'ptbcal');
    calFile = 'SonyPVM2541';
    PTBcal = LoadCalFile(calFile, [], calDir);
    
    %% Force pixel pitch to be the same in horiz and vert directions, because ISETBIO just stores one number
    originalPTBdotsPerMeterX = PTBcal.describe.displayDescription.screenSizePixel(1)/(PTBcal.describe.displayDescription.screenSizeMM(1)/1000);
    PTBcal.describe.displayDescription.screenSizeMM(2) = 1000*PTBcal.describe.displayDescription.screenSizePixel(2)/originalPTBdotsPerMeterX;
    originalPTBdotsPerMeterY = PTBcal.describe.displayDescription.screenSizePixel(2)/(PTBcal.describe.displayDescription.screenSizeMM(2)/1000);

    %% Generate isetbio display object
    extraData = ptb.ExtraCalData;
    extraData.distance = 0.764;
    saveDisplayObject = false;
    isetbioDisplayObject = ptb.GenerateIsetbioDisplayObjectFromPTBCalStruct('SonyPVM2541', PTBcal, extraData, saveDisplayObject);

    %% Now go back into PTB format from the isetbio version
    reconstructedPTBcal  = ptb.GeneratePTCalStructFromIsetbioDisplayObject(isetbioDisplayObject);
    
    %% Comparisons
    % 
    % Display size
    tolerance = 1e-6;
    UnitTest.assertIsZero(abs(PTBcal.describe.displayDescription.screenSizeMM-reconstructedPTBcal.describe.displayDescription.screenSizeMM),'Display size mm',tolerance);
    UnitTest.assertIsZero(abs(PTBcal.describe.displayDescription.screenSizePixel-reconstructedPTBcal.describe.displayDescription.screenSizePixel),'Display size pixel',tolerance);

    % Display resolution (along horizontal dimension)
    %
    % For this test, we'll make the PTB dpi the same for horizontal and
    % vertical, because ISETBIO just takes one number and we can't get two
    % different numbers back.
    isetbioDotsPerMeter = displayGet(isetbioDisplayObject, 'dots per meter');
    reconstructedPTBdotsPerMeterX = reconstructedPTBcal.describe.displayDescription.screenSizePixel(1)/(reconstructedPTBcal.describe.displayDescription.screenSizeMM(1)/1000);
    UnitTest.assertIsZero(abs(originalPTBdotsPerMeterX-reconstructedPTBdotsPerMeterX),'Horizontal resolution',tolerance);
    
    % Display resolution (along vertical dimension)
    isetbioDotsPerMeter = displayGet(isetbioDisplayObject, 'dots per meter');
    reconstructedPTBdotsPerMeterY= reconstructedPTBcal.describe.displayDescription.screenSizePixel(2)/(reconstructedPTBcal.describe.displayDescription.screenSizeMM(2)/1000);
    UnitTest.assertIsZero(abs(originalPTBdotsPerMeterY-reconstructedPTBdotsPerMeterY),'Vertical resolution',tolerance);
    
    % Wavelength sampling
    UnitTest.assertIsZero(abs(PTBcal.describe.S-reconstructedPTBcal.describe.S),'Device gamma table',0);

    % nDevices
    UnitTest.assertIsZero(abs(PTBcal.nDevices-reconstructedPTBcal.nDevices),'Device gamma table',0);

    % Primaries
    UnitTest.assertIsZero(abs(PTBcal.P_device-reconstructedPTBcal.P_device),'Device primary spectra',0);
    
    % Ambient
    UnitTest.assertIsZero(abs(PTBcal.P_ambient-reconstructedPTBcal.P_ambient),'Device ambient spectrum',0);
    
    % Gamma table
    UnitTest.assertIsZero(abs(PTBcal.gammaInput-reconstructedPTBcal.gammaInput),'Device gamma table',0);
    UnitTest.assertIsZero(abs(PTBcal.gammaTable-reconstructedPTBcal.gammaTable),'Device gamma table',0);

    %% Plots
    if (runTimeParams.generatePlots)
    end
    
end