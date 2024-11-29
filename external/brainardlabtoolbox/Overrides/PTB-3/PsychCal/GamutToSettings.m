function settings = GamutToSettings(calOrCalStruct, gamut)
% settings = GamutToSettings(calOrCalStruct, gamut)
%
% Find the best device settings to produce
% the passed linear device coordinates.
% 
% The passed coordinates should be in the range [0,1].
% The returned settings also run from [0,1], but after
% inversion of the device's gamma measurements.
%
% The returned argument values is what you actually should
% get after quantization error.
%
% 9/26/93    dhb   Added calData argument.
% 10/19/93   dhb   Allow gamma table dimensions to exceed device settings.
% 11/11/93   dhb   Update for new calData routines.
% 8/4/96     dhb   Update for stuff bag routines.
% 8/21/97    dhb   Update for structures.
% 4/13/02	 awi   Replaced SettingsToDevice with new name SettingsToPrimary.
% 11/16/06   dhb   Adjust for [0,1] world.  Involves changing what's passed
%                  in and out
% 5/26/12    dhb   Add gammaMode == 2 case.
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


% Retrieve necessary fields from calStruct 
gammaInput   = calStructOBJ.get('gammaInput');
gammaTable   = calStructOBJ.get('gammaTable');
iGammaTable  = calStructOBJ.get('iGammaTable');
gammaMode    = calStructOBJ.get('gammaMode');

% Error checking
if isempty(gammaTable)
	error('No gamma table present in calibration structure');
end
if isempty(gammaMode)
	error('SetGammaMethod has not been called on calibration structure');
end

if gammaMode == 0
	settings = GamutToSettingsSch(gammaInput, gammaTable, gamut);
elseif gammaMode == 1
	if isempty(iGammaTable)
		error('Inverse gamma table not present for gammaMode == 1');
	end
	settings = GamutToSettingsTbl(iGammaTable, gamut);
elseif gammaMode == 2
    settings = GamutToSettingsExhaust(gammaInput, gammaTable, gamut);
else
	error('Requested gamma inversion mode %g is not yet implemented', gammaMode);
end
