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
%    obj            - The scene3D object to render
%    varargin       - (Optional) Other key/value pair arguments
%
% Outputs:
%    ieObject       - The Optical Image object
%    terminalOutput - Terminal output
%

%%
p = inputParser;
p.addRequired('obj',@(x)(isa(x,'sceneEye')));
p.addParameter('scaleIlluminance',true,@islogical);

p.parse(obj,varargin{:});
scaleIlluminance = p.Results.scaleIlluminance;

%% Write out into a pbrt file
objNew = obj.write();
recipe = objNew.recipe; % Update the recipe within the sceneEye object. 

%% Render the pbrt file using docker
%scaleFactor = [];
[ieObject, terminalOutput] = piRender(recipe,'version',3);
        
%% Set OI parameters correctly:
if(~obj.debugMode)
    ieObject = obj.setOI(ieObject,'scaleIlluminance',scaleIlluminance);
end

end
