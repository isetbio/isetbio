function [varargout] = SetSensorColorSpace(calOrCalStruct, T_sensor, S_sensor,quiet)
% Usage (new style):
% SetSensorColorSpace(calOrCalStruct, T_sensor, S_sensor,[quiet])
% where calOrCalStruct is a @CalStruct object.
% To obtain the status of the operation use:
% calOrCalStruct.get('M_linear_device_status')
%
% Usage (old style):
% [cal, M_linear_device_status] = SetSensorColorSpace(calOrCalStruct, T_sensor, S_sensor, [quiet])
% where cal and calOrCalStruct are both old-style cal structures.
%
% Initialize the sensor color space for use in calibration. First input
% argument must be either a @CalStruct object (new style) or a calibration 
% structure (old style). For more info on @CalStruct: 'doc CalStruct'. 
% If that input is neither of these class types, or if the input does not contain 
% the required fields, an error is thrown.
%
% Checks that wavelength sampling is consistent and splines
% if not.
%
% M_linear_device_status indicates status of the operation.
%   == 0: OK
%   == 1: Underdetermined or near zero linear/device conversion matrix
%   == 2: Overdetermined linear/device conversion matrix
%
% quiet flag suppresses error messages, default 1.
%
% 9/13/93     dhb   Wrote it.
% 10/16/93    jms   Added optional calData arg which is freed if passed.
% 10/30/93    dhb   Added nDevices, nBases arguments.
%             dhb   Removed optional freeing.  It is a little too clever.
% 11/11/93    dhb   Re-wrote entirely in Matlab.  The C version is too clever.
% 11/17/93    dhb   Call back through new, less clever, C routines.
%             dhb   Support to store M_ambient_linear with compute vars
% 8/4/96      dhb   Converted to modern format from original DHB scheme.
% 8/21/97     dhb   Converted for structures.
% 3/10/98	  dhb	Store T_linear and S_linear.  Seems like a good idea.
%		            Remove nBasesIn, which is never used.
% 					Change nBasesOut to nPrimaryBases.
% 1/4/00      mpr   Added quiet flag to suppress display of messages I don't care about.
% 2/25/99     dhb   Fix case of SplineCmf.
% 4/5/02      dhb, ly  New calling convention.  Internal naming not updated.
% 4/23/04     dhb   Make quiet the default.
% 3/18/10     dhb   Store T_sensor, S_sensor in fields cal.T_sensor, cal.S_sensor.
%                   These are redundant with old cal.T_linear, cal.S_linear, but
%                   it's possible that deleting those will break something in
%                   some calling program.  So I'm leaving both in for the next
%                   few years.
% 5/08/14     npc   Modifications for accessing calibration data using a @CalStruct object.
%                   The first input argument can be either a @CalStruct object (new style), or a cal structure (old style).
%                   Passing a @CalStruct object is the preferred way because it results in 
%                   (a) less overhead (@CalStruct objects are passed by reference, not by value), and
%                   (b) better control over how the calibration data are accessed.

% Set default value for quiet.
if nargin < 4 || isempty(quiet)
  quiet = 1;
end

% Specify @CalStruct object that will handle all access to the calibration data.
[calStructOBJ, inputArgIsACalStructOBJ] = ObjectToHandleCalOrCalStruct(calOrCalStruct);
if (inputArgIsACalStructOBJ)
    % The input (calOrCalStruct) is a @CalStruct object. Make sure that SetSensorColorSpace is called with no return variables.
    if (nargout > 0)
        error('There should be NO return parameters when calling SetSensorColorSpace with a @CalStruct input. For more info: doc SetSensorColorSpace.');
    end
else
    % The input (calOrCalStruct) is a cal struct. Clear it to avoid  confusion.
    clear 'calOrCalStruct';
end
% From this point onward, all access to the calibration data is accomplised via the calStructOBJ.


% Extract needed colorimetric data
P_device        = calStructOBJ.get('P_device');
T_device        = calStructOBJ.get('T_device');
S               = calStructOBJ.get('S');
nDevices        = calStructOBJ.get('nDevices');
nPrimaryBases   = calStructOBJ.get('nPrimaryBases');

if isempty(P_device) || isempty(T_device) || isempty(S) || ...
        isempty(nDevices) || isempty(nPrimaryBases)
    error('Calibration structure does not contain device colorimetric data');
end

% Extract needed ambient data
P_ambient = calStructOBJ.get('P_ambient');
T_ambient = calStructOBJ.get('T_ambient');
S = calStructOBJ.get('S');
if isempty(P_ambient) || isempty(T_device) || isempty(S)
	error('Calibration structure does not contain ambient data');
end

if CheckWls(S,S_sensor,quiet)
  if ~quiet
    disp('InitCal: Splining T_sensor to match T_device');
  end
  T_sensor = SplineCmf(S_sensor,T_sensor,S);
end

% Compute conversion matrix between device and linear coordinates
M_tmp = M_TToT(T_device,T_sensor);
M_device_linear = M_tmp*P_device;

% Pull out only requested nDevices columns of M_device_linear.
% This is a simple way to define a matrix that maps to device
% space.  More sophisticated methods are possible, which is 
% why we carry the nBasesIn variable around.  
tmp_M_device_linear = M_device_linear(:,1:nDevices);
[m,n] = size(tmp_M_device_linear);

% M_linear_device_status
M_linear_device_status = 0;

if (cond(tmp_M_device_linear) > 1e7) || (m > n)
  M_linear_device_status = 1;
  M_linear_device = pinv(tmp_M_device_linear);
elseif n > m
  M_linear_device_status = 2;
  M_linear_device = pinv(tmp_M_device_linear);
else
  M_linear_device = inv(tmp_M_device_linear);
end

% Convert ambient to linear color space
[M_ambient_linear] = M_TToT(T_ambient,T_sensor);
ambient_linear = M_ambient_linear*P_ambient;

% Update calStructOBJ with computed values
calStructOBJ.set('T_sensor',                T_sensor);           
calStructOBJ.set('S_sensor',                S_sensor);           
calStructOBJ.set('T_linear',                T_sensor);           
calStructOBJ.set('S_linear',                S_sensor);           
calStructOBJ.set('M_device_linear',         M_device_linear);   
calStructOBJ.set('M_linear_device',         M_linear_device);
calStructOBJ.set('M_linear_device_status',  M_linear_device_status);
calStructOBJ.set('M_ambient_linear',        M_ambient_linear);   
calStructOBJ.set('ambient_linear',          ambient_linear);     

if (~inputArgIsACalStructOBJ)
    % Old-style functionality. Return modified cal and M_linear_device_status.
    varargout{1} = calStructOBJ.cal;
    varargout{2} = M_linear_device_status;
    % calStructOBJ is not needed anymore. So clear it from the memory.
    clear 'calStructOBJ'
end
