function d = displayCreate(displayName, varargin)
% Create a display structure
%
%  d = displayCreate(displayFileName, varargin)
%
% Display (d) calibration data are stored in a display structure. These are
% the spectral radiance distribution of its primaries and a gamma function.
%
% Inputs:
%   displayName - Name of a file containing a calibrated display structure.
%                 The supported display files are stored in data/displays.
%                 The files should contain a variable ('d') as display
%                 structure. See displayGet and displaySet for the slots.
%   varargin    - User defined parameter values, should be in name-value
%                 pairs. See displaySet for supported parameters
%
% Example:
%   d = displayCreate;
%   d = displayCreate('lcdExample');
%   wave = 400:5:700; d = displayCreate('lcdExample', 'wave', wave);
%
%   Some displays have psf data, as well, e.g. 
%     d = displayCreate('LCD-Apple');
%
% See Also:
%   sceneFromFile, displayGet, displaySet
%  
% HJ, ISETBIO TEAM, 2015


%% Init Parameters
% Default changed on Nov. 30, 2015.  The original default was far too
% bright for common practice.
if notDefined('displayName'), displayName = 'LCD-Apple'; end

% Identify the object type
d.type = 'display';

% This will change the filename to lower case which can cause problems.
% displayName = ieParamFormat(displayName);

d = displaySet(d, 'name', displayName);

% We can create some displays, or if it is not on the list perhaps it is a
% file name that we load.
switch displayName
    case 'default'
        % See comment about the default above.  We should make it a little
        % closer to sRGB standard chromaticities.
        d = displayDefault(d);
 
    otherwise
        % Read a file with calibrated display data.
        % This can include pixel psf data for some displays.
        if exist(displayName,'file') || exist([displayName,'.mat'],'file') 
            tmp = load(displayName);
            if ~isfield(tmp,'d')
                error('No display struct in the file');
            else  d = tmp.d;
            end
        else error('Unknown display %s.',displayName);
        end

end

if length(varargin) == 1
    warning('ISETBIO: Should set wave as name-value pairs');
    d = displaySet(d, 'wave', varargin{1});
else
    assert(~isodd(length(varargin)), 'varargin should in pairs');
    for ii = 1 : 2 : length(varargin)
        d = displaySet(d, varargin{ii}, varargin{ii+1});
    end
end

% Set the default scene rgb for ISETBIO
d.mainimage = sceneGet(sceneCreate,'rgb image'); % Image for main display window

end

% Create a default display structure
function d = displayDefault(d)
% Create a default display that works well with the imageSPD rendering
% routine.  See vcReadImage for more notes.  Or move those notes here.
%
wave = 400:10:700;
spd = pinv(colorBlockMatrix(length(wave)))/700;   % Makes peak about 100 cd/m2
d = displaySet(d, 'wave', wave);
d = displaySet(d, 'spd', spd);

% Linear gamma function
N = 256; % 8 bit display
g = repmat(linspace(0, 1, N),3,1)';
d = displaySet(d, 'gamma', g);  % From digital value to linear intensity

% Spatial matters
d.dpi = 96;    % Typical display density
d.dist = 0.5;  % Typical viewing distance, 19 inches


end