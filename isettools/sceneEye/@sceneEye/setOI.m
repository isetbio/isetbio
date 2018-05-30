function [ieObject] = setOI(obj,ieObject,varargin)
% Set optical image parameters to match sceneEye parameters. 
%
% Syntax:
%   [success] = write(obj, [varargin])
%
% Description:
%     Given an ieObject (the optical image), we will set the parameters within
%     to match those given in the sceneEye object. This is necessary, since the
%     optical image returned by piRender will have parameters such as f-number,
%     focal length, etc. set to default, but we want these values to reflect
%     the retinal image we have just rendered. 
%
% Inputs:
%    obj            - The scene3D object to render
%    ieObject       - the rendered optical image; typically from piRender.
%
% Outputs:
%   ieObject        - the optical image is returned with the correct
%                     parameters set. 
%

    
ieObject = oiSet(ieObject,'name',sprintf('%s-%s',obj.name,datestr(now,'mmm-dd,HH:MM')));

% Scene distance. We set it to infinity, since it doesn't apply to
% the raytracing.
ieObject = oiSet(ieObject, 'distance', Inf);

% For the optics focal length, we use the distance between the back of
% the lens and the retina. Although it is not exactly the same as the
% focal length for the eye (which also changes with accommodation), it
% is a good approximation for the oiWindow.
ieObject = oiSet(ieObject, 'optics focal length', ...
    obj.retinaDistance * 1e-3);
ieObject = oiSet(ieObject, 'optics fnumber', ...
    obj.retinaDistance / obj.pupilDiameter);
ieObject = oiSet(ieObject, 'fov', obj.fov);

% Clear default optics that do not apply to the iset3d optical
% image. We may want to add these in in the future.
ieObject.optics = opticsSet(ieObject.optics, 'model', 'iset3d');
ieObject.optics = opticsSet(ieObject.optics, 'name', ...
    'PBRT Navarro Eye');
ieObject.optics.OTF = [];

% BW:  We should set the lens density to the value used in
% sceneEye, not just remove it.  Ask TL whether she does anything
% with the lens at all ... if she doesn't, we might apply it here.
% ieObject.optics.lens = [];
ieObject.optics.lens.name = obj.recipe.get('lens file');
ieObject.optics.offaxis = '';
ieObject.optics.vignetting = [];

% Shouldn't we adjust the mean illuminance to some reasonable
% level here?
% disp('myScene.render: Using oiAdjustIlluminance to set mean illuminance to 5 lux.');
% ieObject = oiAdjustIlluminance(ieObject,5);  % 5 lux

end
