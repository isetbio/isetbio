function obj = piWRS(SE,varargin)
% Calls piWRS
%
% Synopsis
%   obj = sceneEye.piWRS(varargin);
%
% Brief
%   Writes, Renders, and Shows the sceneEye recipe
%
% Inputs
%   SE - sceneEye object
%
% Key/val pairs
%    render type      - Usual recipe render type cell array 
%    scaleIlluminance - Scale the returned illuminance with respect to
%                       the pupil area. Default true. 
%    write - Call piWrite first. Default: true 
%            For debugging we sometimes suppress piWrite. 
%    docker wrapper - Specify a docker wrapper.  If not specified, we use
%                         dockerWrapper.humanEyeDocker
%   'name'  - Set the Scene or OI name
%   'gamma' - Set the display gamma for the window
%   'show'  - Default: true
%
% Outputs
%   obj - Returns the ISETBio/Cam scene or oi struct
%
% See also
%   sceneEye.render, sceneEye.write
%

%% Parse
varargin = ieParamFormat(varargin);

p = inputParser;
p.KeepUnmatched = true;
p.addRequired('SE', @(x)(isa(x, 'sceneEye')));
p.addParameter('scaleilluminance', true, @islogical);
p.addParameter('dockerwrapper',[],@(x)(isa(x,'dockerWrapper') || isempty(x)));

p.parse(SE, varargin{:});
scaleIlluminance  = p.Results.scaleilluminance;
thisDockerWrapper = p.Results.dockerwrapper;

thisR = SE.recipe;

% For ISETBio debugging, we sometimes switch the camera to pinhole
if SE.usePinhole
    % We will render a scene through a pinhole camera.  We try to match the
    % fov for the scene with the fov that was set for the eyeballc ase.
    fov = thisR.get('fov');
    cameraSave = thisR.get('camera');
    
    thisR.set('camera',piCameraCreate('pinhole'));
    thisR.set('fov',fov);
end

% Call the master piWRS routine
obj = piWRS(thisR,varargin{:},'dockerwrapper',thisDockerWrapper);

% Deal with special ISETBio pinhole management
if(~SE.usePinhole)
    % If we are not in debug mode with a pinhole, set OI parameters.
    obj = obj.setOI(obj, 'scale illuminance', scaleIlluminance);
    % oiWindow(obj);
else
    % If debugMode, put back the saved camera information.
    thisR.set('camera',cameraSave);
    % sceneWindow(obj);
end

end

%{
varargin = ieParamFormat(varargin);
p = inputParser;

p.addRequired('SE',@(x)(isa(x,'sceneEye')));

p.addParameter('dockerwrapper',[],@(x)(isa(x,'dockerWrapper')));
p.addParameter('name','',@ischar);
p.addParameter('show',true,@islogical);
p.addParameter('gamma',[],@isnumeric);
p.addParameter('renderflag','',@ischar);

p.parse(SE,varargin{:});

thisDocker  = p.Results.dockerwrapper;
g           = p.Results.gamma;
name        = p.Results.name;
show        = p.Results.show;
renderFlag  = p.Results.renderflag;

if isempty(thisDocker)
    thisDocker = dockerWrapper;
    thisDocker.preset('human eye');
end

%% Render

obj = SE.render('dockerwrapper',thisDocker);

%% Edit parameters

switch lower(obj.type)
    case 'scene'

        if ~isempty(name), obj = sceneSet(obj,'name',name); end
        if show
            sceneWindow(obj);
            if ~isempty(g), sceneSet(obj,'gamma',g); end
            if ~isempty(renderFlag), sceneSet(obj,'render flag',renderFlag); end
        end

    case 'opticalimage'

        if ~isempty(name), obj = oiSet(obj,'name',name); end
        if show
            oiWindow(obj);
            if ~isempty(g), oiSet(obj,'gamma',g); end
            if ~isempty(renderFlag), oiSet(obj,'render flag',renderFlag); end
        end

    otherwise
        error('Unknown ISET Object type %s',obj.type)
end

end
%}
