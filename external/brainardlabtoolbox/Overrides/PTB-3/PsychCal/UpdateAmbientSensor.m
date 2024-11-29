function [varargout] = UpdateAmbientSensor(calOrCalStruct,new_ambient_sensor,ADD)
% Usage (new style):
% UpdateAmbientSensor(calOrCalStruct,new_ambient_sensor,ADD)
% where calOrCalStruct is a @CalStruct object.
%
% Usage (old style):
% cal = UpdateAmbientSensor(calOrCalStruct,new_ambient_sensor,[ADD])
% where cal and calOrCalStruct are both old-style cal structures.
%
% Update the ambient light used in the conversions.  The
% value for new_ambient_sensor should be passed in the
% same units as defined by T_sensor in the call to
% SetColorSpace.
%
% If flag ADD is true, passed ambient is added to current
% value.  Otherwise passed value replaces current value.
% ADD is false if not passed.  Use caution when setting ADD
% true -- if the ambient is changing during the experiment
% you typically don't want to keep adding multiple times.
%
% If instead you want to update in the measurement units,
% call UpdateAmbient instead.
%
% 7/7/98    dhb		   Wrote it.
% 4/5/02    dhb, ly  Update for new interface.  Internal names not changed.
% 5/2/02    dhb, kr  Add ADD flag.
% 5/08/14   npc      Modifications for accessing calibration data using a @CalStruct object.
%                    The first input argument can be either a @CalStruct object (new style), or a cal structure (old style).
%                    Passing a @CalStruct object is the preferred way because it results in 
%                    (a) less overhead (@CalStruct objects are passed by reference, not by value), and
%                    (b) better control over how the calibration data are accessed.
%

    % Specify @CalStruct object that will handle all access to the calibration data.
    [calStructOBJ, inputArgIsACalStructOBJ] = ObjectToHandleCalOrCalStruct(calOrCalStruct);
    if (inputArgIsACalStructOBJ)
        % The input (calOrCalStruct) is a @CalStruct object. Make sure that UpdateAmbientSensor is called with no return variables.
        if (nargout > 0)
            error('There should be NO return parameters when calling UpdateAmbientSensor with a @CalStruct input. For more info: doc UpdateAmbientSensor.');
        end
    else
         % The input (calOrCalStruct) is a cal struct. Clear it to avoid  confusion.
        clear 'calOrCalStruct';
    end
    % From this point onward, all access to the calibration data is accomplised via the calStructOBJ.

    % Primitive dimension check
    if (size(new_ambient_sensor,1) ~= size(calStructOBJ.get('ambient_linear'),1) || ...
      size(new_ambient_sensor,2) ~= size(calStructOBJ.get('ambient_linear'),2) )
        error('Old and new ambient specifications are not of same dimension');
    end

    % Set default on optional argument.
    if (nargin < 3 || isempty(ADD))
        ADD = 0;
    end

    % Update the structure
    if (~ADD)
        new_ambient_linear = new_ambient_sensor;
    else
        new_ambient_linear = calStructOBJ.get('ambient_linear') + new_ambient_sensor;
    end

    calStructOBJ.set('ambient_linear', new_ambient_linear);

    if (~inputArgIsACalStructOBJ)
        % Old-style functionality. Return modified cal.
        varargout{1} = calStructOBJ.cal;
        % calStructOBJ is not needed anymore. So clear it from the memory.
        clear 'calStructOBJ'
    end
    
end


