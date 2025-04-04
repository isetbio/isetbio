function PTBcal = GenerateEmptyPTBcalStruct
% Generate an empty PTB calibration structure
%
% Synopsis:
%   PTBcal = ptb.GenerateEmptyPTBcalStruct
%
% Description:
%   Generate an empty PTB calibration structure.
%
% Inputs:
%    None.
%
% Outputs:
%    PTBCal    - The PTB calibration struct.
%
% Optional key/value pairs:
%    None.
%
% See also: ptb.GenerateIsetbioDisplayObjectFromPTBCalStruct,
%           ptb.GeneratePTCalStructFromIsetbioDisplayObject,
%           ptb.GeneratePsychToolboxCalStruct, ptb.GenerateEmptyPTBStruct,
%           generateCustomDisplay.
%

% History:
%  1/16/22 dhb Added comments.

displayDescriptionStruct = struct( ...
     'screenSizeMM', [], ...
  'screenSizePixel', [], ...
      'refreshRate', [], ...
     'bitsPerPixel', [], ...
    'bitsPerSample', [], ...
  'samplesPerPixel', [] ...
  );
    
    gammaStruct = struct(...
            'fitType', 'crtLinear', ...
     'contrastThresh', 0.001, ...
          'exponents', [], ...
     'fitBreakThresh', [], ...
          'useweight', [] ...
    );
    
    describeStruct = struct( ...
                     'S', [], ...
      'blankOtherScreen', [], ...
         'blankSettings', [], ...
            'boxOffsetX', [], ...
            'boxOffsetY', [], ...
               'boxSize', [], ...
       'calibrationType', '', ...
               'caltype', '', ...
               'comment', '', ...
              'computer', '', ...
               'dacsize', [], ...
                  'date', '', ...
    'displayDescription', displayDescriptionStruct, ...
                'driver', '', ...
                 'gamma', gammaStruct, ...
                    'hz', [], ...
         'leaveRoomTime', [], ...
         'meterDistance', [], ...
               'monitor', '', ...
              'nAverage', [], ...
                 'nMeas', [], ...
               'program', '', ...
         'promptforname', [], ...
       'screenSizePixel', [], ...
               'svnInfo', struct(), ...
      'whichBlankScreen', [], ...
        'whichMeterType', [], ...
           'whichScreen', [], ...
                   'who', '', ...
           'yokedmethod', [] ...
        );
    
    
    rawDataStruct = struct(...
       'rawGammaInput', [], ...
       'rawGammaTable', [] ...
    );

    % Generate a calStruct suitable for use with PTB Cal functions
    PTBcal = struct(...
            'describe', describeStruct, ...
    'M_ambient_linear', [], ...
     'M_device_linear', [], ...
     'M_linear_device', [], ...
           'P_ambient', [], ... 
            'P_device', [], ...
           'S_ambient', [], ...
            'S_device', [], ...
            'S_linear', [], ...
            'S_sensor', [], ...
           'T_ambient', [], ...
            'T_device', [], ...
            'T_linear', [], ...
            'T_sensor', [], ...
      'ambient_linear', [], ...
           'basicmeas', struct(), ...
             'bgColor', [], ...
              'bgmeas', struct(), ...
             'fgColor', [], ...
         'gammaFormat', [], ...
          'gammaInput', [], ...
           'gammaMode', [], ...
          'gammaTable', [], ...
         'iGammaTable', [], ...
             'rawdata', rawDataStruct, ...
            'nDevices', [], ...
       'nPrimaryBases', [] ...
     );

end
