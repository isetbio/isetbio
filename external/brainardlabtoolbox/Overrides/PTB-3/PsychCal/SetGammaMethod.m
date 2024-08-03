function [varargout] = SetGammaMethod(calOrCalStruct,gammaMode,precision)
% Usage (new style):
% SetGammaMethod(calOrCalStruct,gammaMode,[precision])
% where calOrCalStruct is a @CalStruct object.
%
% Usage (old style):
% cal = SetGammaMethod(calOrCalStruct,gammaMode,[precision])
% where cal and calOrCalStruct are both old-style cal structures.
%
% Set up the gamma correction mode to be used.  Options
% are:
%   gammaMode == 0 - search table using linear interpolation via interp1.
%   gammaMode == 1 - inverse table lookup.  Fast but less accurate.
%   gammaMode == 2 - exhaustive search
%
% If gammaMode == 1, then you may specify the precision of the
% inverse table.  The default is 1000 levels.
%
% See also RefitCalGamma, CalibrateFitGamma, GamutToSettings
%
% 8/4/96  dhb  Wrote it.
% 8/21/97 dhb  Update for structure.
% 3/12/98 dhb  Change name to SetGammaCorrectMode
% 5/26/12 dhb  Add real exhaustive search mode (2). 
% 5/08/14 npc  Modifications for accessing calibration data using a @CalStruct object.
%              The first input argument can be either a @CalStruct object (new style), or a cal structure (old style).
%              Passing a @CalStruct object is the preferred way because it results in 
%              (a) less overhead (@CalStruct objects are passed by reference, not by value), and
%              (b) better control over how the calibration data are accessed.

    % Specify @CalStruct object that will handle all access to the calibration data.
    [calStructOBJ, inputArgIsACalStructOBJ] = ObjectToHandleCalOrCalStruct(calOrCalStruct);
    if (inputArgIsACalStructOBJ)
        % The input (calOrCalStruct) is a @CalStruct object. Make sure that SetGammaMethod is called with no return variables.
        if (nargout > 0)
            error('There should be NO return parameters when calling SetGammaMethod with a @CalStruct input. For more info: doc SetGammaMethod.');
        end
    else
         % The input (calOrCalStruct) is a cal struct. Clear it to avoid  confusion.
        clear 'calOrCalStruct';
    end
    % From this point onward, all access to the calibration data is accomplised via the calStructOBJ.


    % Check that the needed data is available. 
    gammaTable = calStructOBJ.get('gammaTable'); 
    gammaInput = calStructOBJ.get('gammaInput');

    if isempty(gammaTable)
        error('Calibration structure does not contain gamma data');
    end

    % Do the right thing depending on mode.
    if gammaMode == 0
        calStructOBJ.set('gammaMode', gammaMode);
    elseif gammaMode == 1
        if nargin == 2
            precision = 1000;
        end
        iGammaTable = InvertGammaTable(gammaInput,gammaTable,precision);
        calStructOBJ.set('gammaMode', gammaMode);
        calStructOBJ.set('iGammaTable', iGammaTable);
    elseif gammaMode == 2
        calStructOBJ.set('gammaMode', gammaMode);
    else
      error('Requested gamma inversion mode %g is not yet implemented', gammaMode);
    end

    if (~inputArgIsACalStructOBJ)
        % Old-style functionality. Return modified cal.
        varargout{1} = calStructOBJ.cal;
        % calStructOBJ is not needed anymore. So clear it from the memory.
        clear 'calStructOBJ'
    end
end

	
