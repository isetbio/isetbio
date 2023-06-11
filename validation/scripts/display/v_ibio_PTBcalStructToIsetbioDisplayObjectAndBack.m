function v_ibio_PTBcalStructToIsetbioDisplayObjectAndBack(varargin)
%
% Validate conversion of PTBcalStruct to isetbio display object and back.
%
% This compares the PTB reconstruction after conversion into isetbio and
% back with the originally read PTB cal stucture.
%
% See also 
%    v_IrradianceIsomerizations, v_DisplayColorConversion
%

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
dPTB  = PTBcal.describe.displayDescription.screenSizeMM;
dIBIO = reconstructedPTBcal.describe.displayDescription.screenSizeMM;
assert(dPTB(1) - dIBIO(1) < tolerance && dPTB(2) - dIBIO(2) < tolerance);

% Display resolution (along horizontal dimension)
%
% For this test, we'll make the PTB dpi the same for horizontal and
% vertical, because ISETBIO just takes one number and we can't get two
% different numbers back.
% isetbioDotsPerMeter = displayGet(isetbioDisplayObject, 'dots per meter');
dpmIBIO = reconstructedPTBcal.describe.displayDescription.screenSizePixel(1)/(reconstructedPTBcal.describe.displayDescription.screenSizeMM(1)/1000);
assert(abs(originalPTBdotsPerMeterX-dpmIBIO) < tolerance);

% Display resolution (along vertical dimension)
dpmIBIO = reconstructedPTBcal.describe.displayDescription.screenSizePixel(2)/(reconstructedPTBcal.describe.displayDescription.screenSizeMM(2)/1000);
assert(abs(originalPTBdotsPerMeterY-dpmIBIO) < tolerance);

%     % Wavelength sampling
%     UnitTest.assertIsZero(abs(PTBcal.describe.S-reconstructedPTBcal.describe.S),'Device gamma table',0);
%
%     % nDevices
%     UnitTest.assertIsZero(abs(PTBcal.nDevices-reconstructedPTBcal.nDevices),'Device gamma table',0);
%
%     % Primaries
%     UnitTest.assertIsZero(abs(PTBcal.P_device-reconstructedPTBcal.P_device),'Device primary spectra',0);
%
%     % Ambient
%     UnitTest.assertIsZero(abs(PTBcal.P_ambient-reconstructedPTBcal.P_ambient),'Device ambient spectrum',0);
%
%     % Gamma table
%     UnitTest.assertIsZero(abs(PTBcal.gammaInput-reconstructedPTBcal.gammaInput),'Device gamma table',0);
%     UnitTest.assertIsZero(abs(PTBcal.gammaTable-reconstructedPTBcal.gammaTable),'Device gamma table',0);
end