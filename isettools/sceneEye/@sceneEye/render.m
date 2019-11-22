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

%%
p = inputParser;
p.addRequired('obj', @(x)(isa(x, 'sceneEye')));
p.addParameter('scaleIlluminance', true, @islogical);
p.addParameter('reuse', false, @islogical);

p.parse(obj, varargin{:});
scaleIlluminance = p.Results.scaleIlluminance;
reuse = p.Results.reuse;

%% Write out into a pbrt file
objNew = obj.write();
recipe = objNew.recipe; % Update the recipe within the sceneEye object.

%% Render the pbrt file using docker
%scaleFactor = [];
if reuse
    [ieObject, terminalOutput] = piRender(recipe, 'version', 3, ...
        'reuse', true);
else
    [ieObject, terminalOutput] = piRender(recipe, 'version', 3);
end

%% Set OI parameters correctly:
if(~obj.debugMode)
    ieObject = obj.setOI(ieObject, 'scaleIlluminance', scaleIlluminance);
end

end
