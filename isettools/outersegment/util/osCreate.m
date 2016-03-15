function obj = osCreate(type)
% Create an outer segment object
%
%   os = osCreate([osType = 'linear']);
%
% There are three osTypes (sublcasses of outer segment objects).  These are
%    @osLinear   -  uses linear impulse response to calculate cone current.
%    @osBioPhys  - uses biophysical model by FMR to calculate cone current.
%    @osIdentity - stores raw RGB data, used for stimulus-referred RGC
%    models.
% See sublcass definitions, e.g. @osLinear/osLinear.m, for more details.
% 
% Inputs:  Subclass type ('linear','biophys','identity')
%
% obj:    the outersegment object
%         The obj defaults are set to turn off the noise except for photon
%         noise.  
%
% See also: @osLinear, @osBioPhys and @osIdentity subclasses for more
%           details of the specific implementations.
%
% Examples:
%   os1 = osCreate();
%   os1 = osCreate('biophys');
%   os1 = osCreate('identity');
%
% JRG ISETBIO Team, Copyright, 2015

p = inputParser; 
p.CaseSensitive = false; p.FunctionName = mfilename;

addRequired( p, 'type');

% Parse and put results into structure p.
p.parse(type); params = p.Results;

%% Create the proper object
switch ieParamFormat(params.type)
    case {'linear','oslinear'}
        obj = osLinear();
    case {'biophys','osbiophys','rieke','osrieke'}
        obj = osBioPhys();
    case {'identity','osidentity'}
        obj = osIdentity();
    otherwise
        obj = osLinear();
end

end

