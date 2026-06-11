function obj = osCreate(type)
% Create an outer segment object
%
% Syntax:
%	os = osCreate([osType]);
%
% Description:
%    There are three osTypes (sublcasses of outer segment objects). These
%    subclasses are:
%       @osLinear   - Use linear impulse response to calculate cone current
%       @osBioPhys  - Use FMR's biophysical model to calculate cone current
%       @osIdentity - Store raw RGB data, for stimulus-referred RGC models
%
%    See sublcass definitions, e.g. @osLinear/osLinear.m, for more details.
%
%    Examples are contained in the code. To access, type 'edit osCreate.m'
%    into the Command Window.
%
% Inputs
%    type - (Optional) Subclass type ('linear', 'biophys', 'identity').
%           Default is linear.
%
% Outputs:
%    obj  - The outersegment object. The obj defaults are set to turn off
%           the noise except for photon noise.
%
% Optional key/value pairs:
%    None.
%
% See Also:
%    @osLinear, @osBioPhys and @osIdentity subclasses for more details of
%    the specific implementations.
%

% History:
%    xx/xx/15  JRG  ISETBIO Team, Copyright, 2015
%    02/12/18  jnm  Formatting

% Examples:
%{
    os1 = osCreate();
    os1 = osCreate('biophys');
    os1 = osCreate('identity');
%}

if notDefined('type'), type = 'linear'; end

p = inputParser; 
p.CaseSensitive = false;
p.FunctionName = mfilename;

addRequired( p, 'type');

% Parse and put results into structure p.
p.parse(type);
params = p.Results;

%% Create the proper object
switch ieParamFormat(params.type)
    case {'linear', 'oslinear'}
        obj = osLinear();
    case {'biophys', 'osbiophys', 'rieke', 'osrieke'}
        obj = osBioPhys();
    case {'identity', 'osidentity'}
        obj = osIdentity();
    case {'displayrgb', 'osdisplayrgb'}
        obj = osDisplayRGB();
    otherwise
        obj = osLinear();
end

end
