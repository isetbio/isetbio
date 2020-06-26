function [ieObject] = setOI(obj, ieObject, varargin)
% Set optical image parameters to match sceneEye parameters. 
%
% Syntax:
%   ieObject = write(obj, [varargin])
%
% Description:
%    Given an ieObject (the optical image), we will set the parameters
%    within to match those given in the sceneEye object. This is necessary,
%    since the optical image returned by piRender will have parameters such
%    as f-number, focal length, etc. set to default, but we want these
%    values to reflect the retinal image we have just rendered.
%
% Inputs:
%    obj                   - Object. The scene3D object to render.
%    ieObject              - Object. the rendered optical image; typically
%                            from piRender.
%
% Outputs:
%   ieObject               - Object. The optical image is returned with the
%                            correct parameters set. 
%
% Optional key/value pairs:
%    meanIlluminancePerMM2 - Numeric. The mean illuminance per mm^2.
%                            Default 5.
%    scaleIlluminance      - Boolean. Whether or not to scale the
%                            illuminance. Default true.
%

%% Initialization
p = inputParser;
varargin = ieParamFormat(varargin);

p.addRequired('obj');
p.addRequired('ieObject');

p.addParameter('meanilluminancepermm2', 5, @isnumeric);
p.addParameter('scaleIlluminance', true, @islogical);

p.parse(obj, ieObject, varargin{:});
meanIlluminancePerMM2 = p.Results.meanilluminancepermm2;
scaleIlluminance = p.Results.scaleIlluminance;

%% Initial calculations
ieObject = oiSet(ieObject, 'name', ...
    sprintf('%s-%s', obj.name, datestr(now, 'mmm-dd, HH:MM')));

% Scene distance. Set it to infinity, since it doesn't apply to raytracing.
ieObject = oiSet(ieObject, 'distance', Inf);

% If there is a crop window we're going to have to adjust the FOV
crop_window = obj.recipe.get('cropwindow');
full_size = 2 * tand(obj.fov / 2) * obj.retinaDistance;
image_size = full_size .* crop_window;

% Convert coordinates so center is (0, 0)
image_size = image_size-(full_size / 2);

% Assume image is square
field_angle(1) = atand(image_size(1) / obj.retinaDistance);
field_angle(2) = atand(image_size(2) / obj.retinaDistance);
fov_crop = abs(field_angle(1) - field_angle(2));

% For the optics focal length, we use the distance between the back of
% the lens and the retina. Although it is not exactly the same as the
% focal length for the eye (which also changes with accommodation), it
% is a good approximation for the oiWindow.
ieObject = oiSet(ieObject, 'optics focal length', ...
    obj.retinaDistance * 1e-3);
ieObject = oiSet(ieObject, 'optics fnumber', ...
    obj.retinaDistance / obj.pupilDiameter);
ieObject = oiSet(ieObject, 'fov', fov_crop);

% Clear default optics that do not apply to the iset3d optical
% image. We may want to add these in in the future.
ieObject.optics = opticsSet(ieObject.optics, 'model', 'iset3d');
ieObject.optics = opticsSet(ieObject.optics, 'name', 'PBRT Navarro Eye');
ieObject.optics.OTF = [];
ieObject.optics.lens.name = obj.recipe.get('lens file');
ieObject.optics.offaxis = '';
ieObject.optics.vignetting = [];

%% Apply lens transmittance
% The following code is from oiCalculateIrradiance.m
ieObject = oiSet(ieObject, 'lens density', obj.lensDensity);
irradiance = oiGet(ieObject, 'photons');
wave = oiGet(ieObject, 'wave');

if isfield(ieObject.optics, 'lens')
    transmittance = opticsGet(ieObject.optics, 'transmittance', ...
        'wave', wave);
else
    transmittance = opticsGet(ieObject.optics, 'transmittance', wave);
end

if any(transmittance(:) ~= 1)
    % Do this in a loop to avoid large memory demand
    transmittance = reshape(transmittance, [1 1 length(transmittance)]);
    irradiance = bsxfun(@times, irradiance, transmittance);
end

ieObject = oiSet(ieObject, 'photons', irradiance);

%% Scale the irradiance with pupil size
% This is already done in piDat2ISET.m. However for the human eye the lens
% file is unique, so we are unable to extract the aperture (at the moment)
% in piDat2ISET.m. Until we add that code, let's make the change here where
% the aperture size is known from the sceneEye object. 
if(scaleIlluminance)
    lensArea = pi * (obj.pupilDiameter / 2) ^ 2;
    meanIlluminance = meanIlluminancePerMM2 * lensArea;
    
    ieObject = oiAdjustIlluminance(ieObject, meanIlluminance);
    ieObject.data.illuminance = oiCalculateIlluminance(ieObject);
end

end
