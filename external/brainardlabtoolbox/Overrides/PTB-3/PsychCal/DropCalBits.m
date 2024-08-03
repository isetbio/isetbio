function [varargout] = DropCalBits(calOrCalStruct,whichScreen,forceBits)
% Usage (new style):
% DropCalBits(calOrCalStruct,whichScreen,[forceBit])
% where calOrCalStruct is a @CalStruct object.
%
% Usage (old style):
% cal = DropCalBits(calOrCalStruct,whichScreen,[forceBit])
% where cal and calOrCalStruct are both old-style cal structures.
%
%
% Drops the bitdepth of a calibration file if
% necessary.  Useful for running programs
% transparently on 8 and 10 bit hardware.
%
% If arg forceBits is passed, it is used as
% the current hardware depth.  Otherwise the
% reported DACBits of whichScreen is used.
%
% This code assumes calibration was done at
% equally spaced levels in RGB settings, as is
% the case with our calibration routines.  May
% not generalize, and I haven't worried about
% the roundoff errors.  Certainly OK for basic
% use.
%
% 2/13/05		dhb		Wrote it.
% 5/28/14       npc     Modifications for accessing calibration data using a @CalStruct object.
%                       The first input argument can be either a @CalStruct object (new style), or a cal structure (old style).
%                       Passing a @CalStruct object is the preferred way because it results in 
%                       (a) less overhead (@CalStruct objects are passed by reference, not by value), and
%                       (b) better control over how the calibration data are accessed.


% Get hardware dac level.  Note that the application
% code should use LoadClut, not SetClut, to access
% full bit depth.
if (nargin > 2 && ~isempty(forceBits))
	hardwareBits = forceBits;
else
	hardwareBits = Screen(whichScreen,'Preference','DACBits');
end

% Specify @CalStruct object that will handle all access to the calibration data.
[calStructOBJ, inputArgIsACalStructOBJ] = ObjectToHandleCalOrCalStruct(calOrCalStruct);
if (inputArgIsACalStructOBJ)
    % The input (calOrCalStruct) is a @CalStruct object. Make sure that SetSensorColorSpace is called with no return variables.
    if (nargout > 0)
        error('There should be NO return parameters when calling SetSensorColorSpace with a @CalStruct input. For more info: doc DropCalBits.');
    end
else
    % The input (calOrCalStruct) is a cal struct. Clear it to avoid  confusion.
    clear 'calOrCalStruct';
end
% From this point onward, all access to the calibration data is accomplised via the calStructOBJ.

% Extract needed data
dacsize       = calStructOBJ.get('dacsize');
nMeas         = calStructOBJ.get('nMeas');

% Force calibration down to 8 bits, which is how we plan to use it.
% Simply refit raw data at correct number of input levels.  
if (dacsize > hardwareBits)
	dacsize         = hardwareBits;
	nInputLevels    = 2^dacsize;
	rawGammaInput   = round(linspace(nInputLevels/nMeas,nInputLevels-1, nMeas))';
    % Update calStructOBJ
    calStructOBJ.set('dacsize', dacsize);
    calStructOBJ.set('rawGammaInput', rawGammaInput);
    CalibrateFitGamma(calStructOBJ);
elseif (dacsize < hardwareBits)
	error('Current hardware has greater bit depth than at calibration.');
end

if (~inputArgIsACalStructOBJ)
    % Old-style functionality. Return modified cal.
    varargout{1} = calStructOBJ.cal;
    % calStructOBJ is not needed anymore. So clear it from the memory.
    clear 'calStructOBJ'
end
    
