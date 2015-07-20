function [newVal,newUnitStr,oldUnitStr] = unitConvert(oldVal,valType,oldUnitStr,newUnitStr)
% [newVal,newUnitStr,oldUnitStr] = unitConvert(oldVal,valType,oldUnitStr)
%
% Value types:
%  'length' - default units: 'm'
%           - recognized units: 'nm', 'um', 'mm', 'cm'
%  'wavelength' - default units: 'nm'
%           - recognized units: 'nm', 'um', 'mm', 'cm'
%  'area'   - default units: 'm2'
%           - recognized units: 'nm2', 'um2', 'mm2', 'cm2'
%  'time'   - default units: 'sec'
%           - recognized units: 'usec', 'msec', 'sec'
%  'irradiance' - default units: 'power/[m2_nm]'
%           - recognized units: 'power/[cm2_nm]'
%
% When oldUnitStr is passed as 'defaultunits', units are taken as the
% default.
%
% Convert units following ISETBO conventions.
%
% NOTE: Wavelength is special cased as a length quantity, because we are so
% deeply steeped in thinking about nm that making the default units for
% wavelength just seems a bridge too far.
%
% NOTE: This routine does not convert units of power from energy to quantal
% units or vice-versa, because that conversion requires knowing the
% wavelength of the light being converted.
%
% 7/20/15  dhb  Started to write this.

% Switch on value type
switch (valueType)
    case 'length'
        % Length
        
        % Factor from old to default
        switch oldUnitStr
            case {'defaultStr', 'm'}
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
                error('Bad units %s passed for type %s',oldUnitStr,valueType);
        end
        
         % Factor from old to default
        switch newUnitStr
            case {'defaultStr', 'm'}
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
                error('Bad units %s passed for type %s',oldUnitStr,valueType);
        end
        
        % Convert
        newVal = oldConverstionFactor*newConversionFactor*oldVal;
    case 'wavelength'
    case 'area'
    case 'time'
    case 'irradiance'
    otherwise
        error('Unknown value type passed');
end


