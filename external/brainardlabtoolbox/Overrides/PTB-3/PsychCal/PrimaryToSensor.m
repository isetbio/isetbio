function [sensor] = PrimaryToSensor(calOrCalStruct,primary)
% [sensor] = PrimaryToSensor(calOrCalStruct,primary)
%
% Convert from primary coordinates to sensor color 
% space coordinates.  The ambient lighting is added to
% the color space coordinates of the device.
%
% This depends on the standard calibration globals.

% 9/26/93    dhb   Added calData argument.
% 10/19/93   dhb   Allow device characterization dimensions to exceed
%                  device settings dimensions.
% 11/11/93   dhb   Update for new calData routines.
% 8/4/96     dhb   Update for new stuff bag routines.
% 8/21/97	 dhb   Convert for structures.
% 4/5/02     dhb   New naming interface.  Internal naming not changed.
% 7/25/10    dhb   Use bsxfun to make it a little faster.
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


% Get necessary calibration data
M_device_linear = calStructOBJ.get('M_device_linear'); 
ambient_linear  = calStructOBJ.get('ambient_linear'); 

if (isempty(M_device_linear) || isempty(ambient_linear))
	error('SetSensorColorSpace has not been called on calibration structure');
end

% Get size
[m,n] = size(primary);

% Color space conversion
[mm,nm] = size(M_device_linear);
if (m > nm)
  error ('Incorrect dimensions for M_device_linear');
end
linear = M_device_linear(:,1:m)*primary;

% Ambient corrections
[ma,na] = size(ambient_linear);
[m,n] = size(linear);
if (ma ~= m || na ~= 1)
  error('Incorrect dimensions for ambient');
end
sensor = bsxfun(@plus,linear,ambient_linear);

