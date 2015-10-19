function obj = osCreate(varargin)
% Create an outer segment object
%
%   os = osCreate([osType = 'linear'],[sensor = null]);
%
% There are three osTypes (sublcasses of outer segment objects).  These are
%    @osLinear   -  Comment me.
%    @osBioPhys  -
%    @osIdentity - 
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

%% Set up the defaults

% We want to use the parser structure if possible (BW)
osType = 'linear'; sensor = [];      % Default values
if ~isempty(varargin),   osType = ieParamFormat(varargin{1}); end
if length(varargin) > 1, sensor = varargin{2}; end

%% Create the proper object
switch osType
    case 'linear'
        obj = osLinear(sensor);
    case {'biophys','rieke'}
        obj = osBioPhys(sensor);
    case 'identity'
        obj = osIdentity(sensor);
    otherwise
        obj = osLinear(sensor);
end

end

