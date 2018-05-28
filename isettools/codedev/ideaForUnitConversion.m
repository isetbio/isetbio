function [newVal, newUnitStr, oldUnitStr] = ideaForUnitConversion(...
    oldVal, valType, oldUnitStr, newUnitStr)
% An idea about how to encapsulate unit conversion.
%
% Syntax:
%   [newVal, newUnitStr, oldUnitStr] = ideaForUnitConversion(...
%       oldVal, valType, oldUnitStr)
%
% Description:
%    This function contains an idea about how to encapsulate unit
%    conversions following ISETBIO conventions.
%
%    When oldUnitStr is passed to the function as 'default' or
%    'defaultStr', the units are taken as the default. Note that converted
%    units follow the existing ISETBIO conventions.
%
% Inputs:
%    oldVal     - Numeric. The original value, in the units specified by
%                 the variable oldUnitStr.
%    valType    - String. A string indicating the type of measurement. The
%                 types are as follows:
%             'length'     - default units: 'm'
%                          - recognized units: 'nm', 'um', 'mm', 'cm'
%             'wavelength' - default units: 'nm'
%                          - recognized units: 'nm', 'um', 'mm', 'cm'
%             'area'       - default units: 'm2'
%                          - recognized units: 'nm2', 'um2', 'mm2', 'cm2'
%             'time'       - default units: 'sec'
%                          - recognized units: 'usec', 'msec', 'sec'
%             'irradiance' - default units: 'power/[m2_nm]'
%                          - recognized units: 'power/[cm2_nm]'
%    oldUnitStr - String. The string specifying the units for the original
%                 value, stored in the variable oldVal.
%    newUnitStr - String. The string specifying the units you wish to
%                 switch to.
%
% Outputs:
%    newVal     - Numeric. The calculated new value, in the units specified
%                 in the variable newUnitStr.
%    newUnitStr - String. The string containing the units of newVal.
%    oldUnitStr - String. The string containing the units associated with
%                 the variable oldVal.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * [Note: XXX - Wavelength is special cased as a length quantity,
%      because we are so deeply steeped in thinking about nm that making
%      the default units for wavelength just seems a bridge too far.]
%    * [Note: XXX - This routine does not convert units of power from
%      energy to quantal units or vice-versa, because that conversion
%      requires knowing the wavelength of the light being converted.]
%

% History:
%    07/20/15  dhb  Started to write this.
%    04/17/18  jnm  Formatting

% Switch on value type
switch valType
    case 'length'
        % Factor from old to default.
        switch oldUnitStr
            case {'default', 'defaultStr', 'm'}
                oldConversionFactor = 1;
            case 'nm'
                oldConversionFactor = 1e-9;
            case 'um'
                oldConversionFactor = 1e-6;
            case 'mm'
                oldConversionFactor = 1e-3;
            case 'cm'
                oldConversionFactor = 1e-2;
            otherwise
                error('Bad units %s passed for type %s', ...
                    oldUnitStr, valueType);
        end

        % Factor from default to new.
        switch newUnitStr
            case {'default', 'defaultStr', 'm'}
                newConversionFactor = 1;
            case 'nm'
                newConversionFactor = 1e9;
            case 'um'
                newConversionFactor = 1e6;
            case 'mm'
                newConversionFactor = 1e3;
            case 'cm'
                newConversionFactor = 1e2;
            otherwise
                error('Bad units %s passed for type %s', ...
                    oldUnitStr, valueType);
        end

    case 'wavelength'
        % Factor from old to default.
        switch oldUnitStr
            case {'default', 'defaultStr', 'm'}
                oldConversionFactor = 1;
            case 'nm'
                oldConversionFactor = 1e-9;
            case 'um'
                oldConversionFactor = 1e-6;
            case 'mm'
                oldConversionFactor = 1e-3;
            case 'cm'
                oldConversionFactor = 1e-2;
            otherwise
                error('Bad units %s passed for type %s', ...
                    oldUnitStr, valueType);
        end

        % Factor from default to new.
        switch newUnitStr
            case {'default', 'defaultStr', 'm'}
                newConversionFactor = 1;
            case 'nm'
                newConversionFactor = 1e9;
            case 'um'
                newConversionFactor = 1e6;
            case 'mm'
                newConversionFactor = 1e3;
            case 'cm'
                newConversionFactor = 1e2;
            otherwise
                error('Bad units %s passed for type %s', ...
                    oldUnitStr, valueType);
        end
    case 'area'
        % Factor from old to default.
        switch oldUnitStr
            case {'default', 'defaultStr', 'm2'}
                oldConversionFactor = 1;
            case 'nm2'
                oldConversionFactor = 1e-18;
            case 'um2'
                oldConversionFactor = 1e-12;
            case 'mm2'
                oldConversionFactor = 1e-6;
            case 'cm2'
                oldConversionFactor = 1e-4;
            otherwise
                error('Bad units %s passed for type %s', ...
                    oldUnitStr, valueType);
        end

        % Factor from default to new.
        switch newUnitStr
            case {'default', 'defaultStr', 'm2'}
                newConversionFactor = 1;
            case 'nm2'
                newConversionFactor = 1e18;
            case 'um2'
                newConversionFactor = 1e12;
            case 'mm2'
                newConversionFactor = 1e6;
            case 'cm2'
                newConversionFactor = 1e4;
            otherwise
                error('Bad units %s passed for type %s', ...
                    oldUnitStr, valueType);
        end
    case 'time'
        % Factor from old to default.
        switch oldUnitStr
            case {'default', 'defaultStr', 'sec'}
                oldConversionFactor = 1;
            case 'usec'
                oldConversionFactor = 1e-6;
            case 'msec'
                oldConversionFactor = 1e-3;
            otherwise
                error('Bad units %s passed for type %s', ...
                    oldUnitStr, valueType);
        end

        % Factor from default to new.
        switch newUnitStr
            case {'default', 'defaultStr', 'sec'}
                newConversionFactor = 1;
            case 'usec'
                newConversionFactor = 1e6;
            case 'msec'
                newConversionFactor = 1e3;
            otherwise
                error('Bad units %s passed for type %s', ...
                    oldUnitStr, valueType);
        end
    case 'irradiance'
        % Factor from old to default.
        switch oldUnitStr
            case {'default', 'defaultStr', 'power/[m2_nm]'}
                oldConversionFactor = 1;
            case 'power/[cm2_nm]'
                oldConversionFactor = 1e-4;
            otherwise
                error('Bad units %s passed for type %s', ...
                    oldUnitStr, valueType);
        end

        % Factor from default to new.
        switch newUnitStr
            case {'default', 'defaultStr', 'power/[m2_nm]'}
                newConversionFactor = 1;
            case 'power/[cm2_nm]'
                newConversionFactor = 1e4;
            otherwise
                error('Bad units %s passed for type %s', ...
                    oldUnitStr, valueType);
        end
    otherwise
        error('Unknown value type passed');
end

% Convert
newVal = oldConversionFactor * newConversionFactor * oldVal;
