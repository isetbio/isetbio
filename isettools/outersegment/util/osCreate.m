function obj = osCreate(varargin)
% Generate an @osLinear, @osBioPhys or @osIdentity object
%
%   os = osCreate([osType = 'linear'],[sensor = null]);
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
% There are three osTypes (sublcasses of outer segment objects).  These are
%    @osLinear   -  uses linear impulse response to calculate cone current.
%    @osBioPhys  - uses biophysical model by FMR to calculate cone current.
%    @osIdentity - stores raw RGB data, used for stimulus-referred RGC
%    models.
% See sublcass definitions, e.g. @osLinear/osLinear.m, for more details.
% 
% Examples:
%   os1 = osCreate();
%   os1 = osCreate('biophys');
%   os1 = osCreate('identity');
%
% JRG ISETBIO Team, Copyright, 2015

%% Set up the defaults

% We want to use the parser structure if possible (BW)
osType = 'linear';   % Default values
if ~isempty(varargin),   osType = ieParamFormat(varargin{1}); end %#ok<*NASGU>

if nargin == 0
    obj = osLinear();
elseif nargin == 1
    if strcmpi(varargin{1},'linear');
        obj = osLinear();
    elseif (strcmpi(varargin{1},'biophys') || strcmpi(varargin{1},'rieke'));
        obj = osBioPhys();
    elseif strcmpi(varargin{1},'identity');
        obj = osIdentity();
    else
        obj = osLinear();
    end
elseif nargin == 3
    if strcmpi(varargin{1},'linear');
        obj = osLinear(varargin{2},varargin{3});
    elseif (strcmpi(varargin{1},'biophys') || strcmpi(varargin{1},'rieke'));
        obj = osBioPhys(varargin{2},varargin{3});
    elseif strcmpi(varargin{1},'identity');
        obj = osIdentity(varargin{2},varargin{3});
    else
        obj = osLinear(varargin{2},varargin{3});
    end
else
    obj = osLinear();
end

end

