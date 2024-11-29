function [varargout] = UpdateAmbient(calOrCalStruct,newP_ambient,ADD)
% Usage (new style):
% UpdateAmbient(calOrCalStruct,newP_ambient,[ADD])
% where calOrCalStruct is a @CalStruct object.
%
% Usage (old style):
% cal = UpdateAmbient(calOrCalStruct,newP_ambient,[ADD])
% where cal and calOrCalStruct are both old-style cal structures.
% 
%
% Update the ambient light used in the conversions.  The
% value of P_ambient in the structure is replaced with the
% passed value and the computed quantities that depend on
% it are updated.
%
% The ambient must be specified in the same measurement units as it
% was in the cal at the initial call to SetColorSpace.  If different
% units are desired, all ambient fields in the structure must be updated
% and SetColorSpace called again.  I have never wanted to do this,
% so I haven't written a separate routine.
%
% It is sometimes useful, however, to update the ambient in the
% linear color space defined by the call to SetColorSpace.  To
% do this, use UpdateAmbientLinear.
%
% If flag ADD is true, passed ambient is added to current
% value.  Otherwise passed value replaces current value.
% ADD is false if not passed.  Use caution when setting ADD
% true -- if the ambient is changing during the experiment
% you typically don't want to keep adding multiple times.
%
% 11/17/93  dhb     Wrote it.
% 8/4/96    dhb     Updated for modern scheme.
% 8/21/97   dhb     Update for structures.
% 3/2/98    dhb		Fix bug in checks introduce 8/21/97, pointed out by dgp.
% 3/10/98	dhb		Change T_ to P_.
% 10/26/99  dhb,mdr Fix bug in checks. There was also a variable name
%					glitch.  I don't think this could have worked the way
%					it was.  Perhaps no one calls it.
% 5/2/02    dhb, kr Add ADD flag.
% 5/28/14   npc     Modifications for accessing calibration data using a @CalStruct object.
%                   The first input argument can be either a @CalStruct object (new style), or a cal structure (old style).
%                   Passing a @CalStruct object is the preferred way because it results in 
%                   (a) less overhead (@CalStruct objects are passed by reference, not by value), and
%                   (b) better control over how the calibration data are accessed.

% Set default on optional argument.
if (nargin < 3 || isempty(ADD))
	ADD = 0;
end

% Specify @CalStruct object that will handle all access to the calibration data.
[calStructOBJ, inputArgIsACalStructOBJ] = ObjectToHandleCalOrCalStruct(calOrCalStruct);
if (inputArgIsACalStructOBJ)
    % The input (calOrCalStruct) is a @CalStruct object. Make sure that UpdateAmbient is called with no return variables.
    if (nargout > 0)
        error('There should be NO return parameters when calling UpdateAmbient with a @CalStruct input. For more info: doc UpdateAmbient.');
    end
else
    % The input (calOrCalStruct) is a cal struct. Clear it to avoid  confusion.
    clear 'calOrCalStruct';
end
% From this point onward, all access to the calibration data is accomplised via the calStructOBJ.


% Extract needed colorimetric data
S                = calStructOBJ.get('S');
oldP_ambient     = calStructOBJ.get('P_ambient');
T_ambient        = calStructOBJ.get('T_ambient');
M_ambient_linear = calStructOBJ.get('M_ambient_linear');

% Check that passed data are compatible
if (isempty(oldP_ambient) || isempty(T_ambient) || isempty(S))
	error('Calibration structure does not contain ambient data');
end

if (isempty(M_ambient_linear))
	error('SetColorSpace has not been called on this calibration structure');
end

[nOld,mOld] = size(oldP_ambient);
[nNew,mNew] = size(newP_ambient);
if (nOld ~= nNew || mOld ~= mNew)
	error('Old and new ambient specifications are not in same units');
end

% Update
if (~ADD)
	P_ambient = newP_ambient;
else
	P_ambient = oldP_ambient + newP_ambient;
end

% Update calStructOBJ with computed values
calStructOBJ.set('ambient_linear',  M_ambient_linear*P_ambient);
calStructOBJ.set('P_ambient', P_ambient);

if (~inputArgIsACalStructOBJ)
    % Old-style functionality.
    varargout{1} = calStructOBJ.cal;
    % calStructOBJ is not needed anymore. So clear it from the memory.
    clear 'calStructOBJ'
end



