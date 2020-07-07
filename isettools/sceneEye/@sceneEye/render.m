function [ieObject, terminalOutput] = render(obj, varargin)
% Render a scene3D object and return an optical image.
%
% Syntax:
%   [ieObject, terminalOutput, outputFile] = render(obj, [varargin])
%
% Description:
%	 Given a scene3D object, we have all the information we need to
%    construct a PBRT file and render it. Therefore, this function does the
%    following:
%       1. Write out a new PBRT ([renderName].pbrt) in working directory
%       2. Render using docker container
%       3. Load the output into an ISETBIO optical image, filling in the
%          right parameters with the scene information
%       4. Return the OI
%
% Inputs:
%    obj              - Object. The scene3D object to render.
%
% Outputs:
%    ieObject         - Object. The Optical Image object.
%    terminalOutput   - String. Terminal output.
%
% Optional key/value pairs:
%    scaleIlluminance - Boolean. Whether or not to calculate the scale
%                       illuminance of the scene.
%    reuse            - Boolean. Whether or not to reuse existing
%                       renderings of the same calculation. (Warning: This
%                       means changes to the parameters will not be
%                       displayed in the rendered image.)
%
% See also
%    piRender

%% Programming
%  Why does this return a scene sometimes and an oi sometimes?
%  

%%
p = inputParser;
p.addRequired('obj', @(x)(isa(x, 'sceneEye')));
p.addParameter('scaleilluminance', true, @islogical);
p.addParameter('reuse', false, @islogical);

p.parse(obj, varargin{:});
scaleIlluminance = p.Results.scaleilluminance;
reuse = p.Results.reuse;

thisR = obj.recipe;

if obj.debugMode
    % We will render a scene through a pinhole camera
    cameraSave = thisR.get('camera');
    thisR.set('camera',piCameraCreate('pinhole'));
    thisR.set('fov',obj.fov);
end

%% Write out into a pbrt file

% Can this just be piWrite(thisR)?  Or does write() do a lot of stuff?

objNew = obj.write();
thisR = objNew.recipe; % Update the recipe within the sceneEye object.

%% Render the pbrt file using docker
%scaleFactor = [];
if reuse
    [ieObject, terminalOutput] = piRender(thisR, 'reuse', true);
else
    [ieObject, terminalOutput] = piRender(thisR);
end

%% If we are not in debug mode, set OI parameters.  
% I guess if we are in debug mode, we return a scene.
if(~obj.debugMode)
    ieObject = obj.setOI(ieObject, 'scale illuminance', scaleIlluminance);
else
    thisR.set('camera',cameraSave);
    obj.recipe = thisR;
end

% oiWindow(ieObject);


end
