function [ieObject, terminalOutput] = render(obj, varargin)
% Render a sceneEye object 
%
% Syntax:
%   [ieObject, terminalOutput] = render(obj, [varargin])
%
% Description:
%	A sceneEye object wraps an iset3d recipe that contains the information 
%   needed to (a) construct a PBRT file, and (b) render it through one of
%   the human physiological optics models.
%
%   This returns a scene when in debugMode or when the camera model is
%   'pinhole' or equivalently 'perspective.
%
%
% Inputs:
%    obj       - Object. The scene3D object to render.  This object
%                has a slot for an iset3d render recipe.
%
% Optional key/value pairs:
%    render type - One of these types of renders.
%      {'radiance','depth','both','all', 
%       'coordinates','material','mesh', 'illuminant','illuminantonly'};  
%         Default is 'both', meaning radiance and depth.
%
%    scaleIlluminance - Boolean. Whether or not to scale the oi
%                       illuminance.
%    write - Typically, we call piWrite() to make sure the pbrt file is
%            updated.  But for debugging we sometimes suppress the piWrite.
%
%
% Outputs:
%    ieObject         - Object. The Optical Image object.
%    terminalOutput   - String. Terminal output.
%
% Description:
%
%   This method renders an iset3d image using physiological optics model.
%   The actions are
%
%    1. Writes out a new PBRT ([renderName].pbrt) in working directory
%    2. Renders using the PBRT spectral docker container
%    3. Loads the output into an ISETBio optical image (unless in
%    debugMode, in which case it is a scene), filling in the parameters
%    with the ISETBio information from the rendering recipe
%
% It returns an oi when there is a lens specified (omni, realisticEye), but
% if you turn on the debugMode it renders a scene through a pinhole.
%
% Dependencies
%   iset3d, ISEBio
%
% See also
%   recipe, piWrite, piRender


%% Parse
varargin = ieParamFormat(varargin);

p = inputParser;
p.addRequired('obj', @(x)(isa(x, 'sceneEye')));
p.addParameter('scaleilluminance', true, @islogical);
p.addParameter('dockerwrapper',[],@(x)(isa(x,'dockerWrapper')));
p.addParameter('write',true,@islogical);

% Some day, check that the cell array has one of these types.
% rTypes = {'radiance','depth','both','all','coordinates','material','mesh', 'illuminant','illuminantonly'};
p.addParameter('rendertype',{'radiance','depth'},@iscell);

p.parse(obj, varargin{:});
renderType        = p.Results.rendertype;
scaleIlluminance  = p.Results.scaleilluminance;
thisDockerWrapper = p.Results.dockerwrapper;

%% Get the render recipe

thisR = obj.recipe;

% For debugging, we sometimes switch the camera to pinhole
if obj.usePinhole
    % We will render a scene through a pinhole camera.  We try to match the
    % fov for the scene with the fov that was set for the eyeballc ase.
    fov = thisR.get('fov');
    cameraSave = thisR.get('camera');
    
    thisR.set('camera',piCameraCreate('pinhole'));
    thisR.set('fov',fov);
end

%% Write out into a pbrt file

if p.Results.write
    % For debugging, we sometimes just render.
    piWrite(thisR);
end

%% Render the pbrt file using docker

% We need a dockerWrapper that works for the human eye when we call
% this.  That should either be set up in the default by setprefs or by
% passing in a specific dockerWrapper.
if isempty(thisDockerWrapper)
    % We decided not to recreate the docker wrapper because we are concerned
    % that the creation resets the docker image on the remote machine
    % and slows things down. If that's false, then it would be fine
    % to recreate thisDockerWrapper as default dockerWrapper above
    % when p.addParameter is called.
    [ieObject, terminalOutput] = piRender(thisR,'render type',renderType);
else
    [ieObject, terminalOutput] = piRender(thisR,'render type',renderType,'ourdocker',thisDockerWrapper);
end

%% Fix up the returned object

if(~obj.usePinhole)
    % If we are not in debug mode with a pinhole, set OI parameters.
    ieObject = obj.setOI(ieObject, 'scale illuminance', scaleIlluminance);
    % oiWindow(ieObject);
else
    % If debugMode, put back the saved camera information.
    thisR.set('camera',cameraSave);
    % sceneWindow(ieObject);
end

end
