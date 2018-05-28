function sFactor = ieUnitScaleFactor(unitName)
% Return scale factor that converts from meters or seconds to other scales
%
% Syntax:
%   sFactor = ieUnitScaleFactor(unitName)
%
% Description:
%    The routine returns the scale factor that converts from meters, or
%    seconds, to other scales. This routine is used in various sceneGet/Set
%    and oiGet/Set operations and oiSpatialSupport. By using this routine,
%    we can specify the units for various returned quantities.
%
%    Examples in code
%
% Inputs:
%    unitName - the unit you wish to convert to. Options include:
%       From meters to:
%           nm     - nanometer(s)
%           um     - micrometer, micron(s)
%           mm     - millimeter(s)
%           cm     - centimeter(s)
%           m      - meter(s)
%           km     - kilometer(s)
%           in     - inch(es)
%           ft     - foot/feet
%       From seconds to:
%           sec    - sec(ond)
%           ms     - millisecond
%           us     - microsecond
%       From radians to:
%           deg    - degrees
%           arcmin - arcmin
%           arcsec - arcsec
%
% Outputs:
%    sFactor  - the requested scale factor
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    11/21/17  jnm  Formatting
%    01/17/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    nMeters = 1;
    nInches = 39.3701;

    % Get number of inches from number of meters
    % Inches is fun for an example.  But never use them.
    nMeters * ieUnitScaleFactor('in')

    % Get number of meters from number of inches
    nInches/ieUnitScaleFactor('in')
%}
%{
    ieUnitScaleFactor('ms')
    ieUnitScaleFactor('deg')
%}

if notDefined('unitName'), error('Unit name must be defined.'); end
unitName = ieParamFormat(unitName);

switch unitName
    % Convert space
    case {'nm', 'nanometer', 'nanometers'}
        sFactor = 1e9;
    case {'micron', 'micrometer', 'um', 'microns'}
        sFactor = 1e6;
    case {'mm', 'millimeter', 'millimeters'}
        sFactor = 1e3;
    case {'cm', 'centimeter', 'centimeters'}
        sFactor = 1e2;
    case {'m', 'meter', 'meters'}
        sFactor = 1;
    case {'km', 'kilometer', 'kilometers'}
        sFactor = 1e-3;
    % Convert meter to English 
    case {'inches', 'inch', 'in'}
        sFactor = 39.37007874;   % inches/meter
    case {'foot', 'feet', 'ft'}
        sFactor = 3.280839895;   % feet/meter
        
    % Convert seconds to other unit 
    case {'s', 'second', 'sec'}
        sFactor = 1;
    case {'ms', 'millisecond'}
        sFactor = 1e3;
    case {'us', 'microsecond'}
        sFactor = 1e6;
        
    % Convert radians to other units
    case {'degrees', 'deg'}
        sFactor = 180 / pi;
    case {'arcmin', 'minutes', 'min'}
        sFactor = (180 / pi) * 60;
    case {'arcsec'}
        sFactor = (180 / pi) * 60 * 60;
        
    otherwise
        error('Unknown spatial unit specification');
end

return;