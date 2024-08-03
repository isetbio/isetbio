function [settings, badIndex] = PrimaryToSettings(calOrCalStruct, primary)
% [settings, badIndex] = PrimaryToSettings(calOrCalStruct, primary)
%
% Convert from primary color space coordinates to device
% setting coordinates.
%
% This depends on the standard calibration globals.

% 9/26/93    dhb   Added calData argument, badIndex return.
% 4/5/02     dhb, ly  New calling interface.
% 5/08/14    npc   Modifications for accessing calibration data using a @CalStruct object.
%                  The first input argument can be either a @CalStruct object (new style), or a cal structure (old style).
%                  Passing a @CalStruct object is the preferred way because it results in 
%                  (a) less overhead (@CalStruct objects are passed by reference, not by value), and
%                  (b) better control over how the calibration data are accessed.

% Specify @CalStruct object that will handle all access to the calibration data.
[calStructOBJ, inputArgIsACalStructOBJ] = ObjectToHandleCalOrCalStruct(calOrCalStruct);
if (~inputArgIsACalStructOBJ)
     % The input (calOrCalStruct) is a cal struct. Clear it to avoid  confusion.
    clear 'calOrCalStruct';
end
% From this point onward, all access to the calibration data is accomplised via the calStructOBJ.


[gamut, badIndex] = PrimaryToGamut(calStructOBJ, primary);
settings = GamutToSettings(calStructOBJ, gamut);
