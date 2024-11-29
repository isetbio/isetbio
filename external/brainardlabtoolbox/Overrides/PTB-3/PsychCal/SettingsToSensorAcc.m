function [sensor,primaryE] = SettingsToSensorAcc(calOrCalStruct,settings)
% [sensor,primaryE] = SettingsToSensorAcc(calOrCalStruct,settings)
%
% Convert from device setting coordinates to
% sensor color space coordinates.  Uses full
% basis function measurements in doing
% conversions so that it can compensate for
% device primary spectral shifts.

% 11/12/93  dhb   Wrote it.
% 11/15/93  dhb   Added deviceE output.
% 8/4/96    dhb   Update for stuff bag routines.
% 8/21/97	dhb	  Update for structures.
% 3/10/98   dhb	  Change nBasesOut to nPrimaryBases.
% 4/5/02    dhb, ly  Update for new calling interface.
% 11/22/09  dhb   Check basis dimension and do the simple fast thing if it is 1.
%                 This will speed things up when there is no point in trying the
%                 fancier algorithm.
% 5/08/14   npc   Modifications for accessing calibration data using a @CalStruct object.
%                 The first input argument can be either a @CalStruct object (new style), or a cal structure (old style).
%                 Passing a @CalStruct object is the preferred way because it results in 
%                 (a) less overhead (@CalStruct objects are passed by reference, not by value), and
%                 (b) better control over how the calibration data are accessed.
   

% Specify @CalStruct object that will handle all access to the calibration data.
[calStructOBJ, inputArgIsACalStructOBJ] = ObjectToHandleCalOrCalStruct(calOrCalStruct);
if (~inputArgIsACalStructOBJ)
    % The input (calOrCalStruct) is a cal struct. Clear it to avoid  confusion.
    clear 'calOrCalStruct';
end
% From this point onward, all access to the calibration data is accomplised via the calStructOBJ.

% Get necessary calibration data
nPrimaryBases = calStructOBJ.get('nPrimaryBases');

if (isempty(nPrimaryBases))
    error('No nPrimaryBases field present in calibration structure');
end

if (nPrimaryBases == 1)
    sensor = SettingsToSensor(calStructOBJ,settings);
    primaryE = [];
else
    settingsE = ExpandSettings(settings,nPrimaryBases);
    primaryE  = SettingsToPrimary(calStructOBJ,settingsE);
    sensor    = PrimaryToSensor(calStructOBJ,primaryE);
end
