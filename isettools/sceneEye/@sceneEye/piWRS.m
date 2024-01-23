function obj = piWRS(SE,varargin)
% Calls the main piWRS function, but accounting for sceneEye
% parameters
%
% Synopsis
%   obj = sceneEye.piWRS(varargin);
%
% Brief
%   Writes, Renders, and Shows the sceneEye recipe.  The typical piWRS
%   parameters are passed through. The special case of accounting for
%   the pinhole testing case and scaling the illuminance case are
%   managed here.  But everything else is passed to the piWRS function
%   in ISET3d-v4.
%
%   The Docker implementation does NOT include multiplication by the
%   lens transmission.  The lens density is stored in the slot
%   sceneEye.lensDensity.  Upon return of the optical image, we create
%   the lens transmission for that lens density and apply it to the
%   OI.
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
% thisDockerWrapper = p.Results.dockerwrapper;

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

% We render but do not show at this point.   We need to apply the
% oi settings before showing.  The varargin can include 'docker
% wrapper' as a key/val pair.
obj = piRender(thisR,varargin{:});

% Deal with special ISETBio pinhole management
if(~SE.usePinhole)
    % If we are not in pinhole (debug) mode, set the OI parameters.
    % This includes applying the lens transmittance.
    obj = SE.setOI(obj, 'scale illuminance', scaleIlluminance);    
else
    % If in pinhole (debug) mode, copy back the saved camera
    % information that was stored above.
    thisR.set('camera',cameraSave);    
end

% Ready to show.
switch obj.type
    case 'opticalimage'
        oiWindow(obj);
    case 'scene'
        sceneWindow(obj);
end

end

