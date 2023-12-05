function [ieObject] = setOI(SE, ieObject, varargin)
% Set optical image parameters to match sceneEye parameters. 
%
% Syntax:
%   ieObject = write(obj, [varargin])
%
% Description:
%    Given an ieObject (the optical image), we set the parameters
%    within to match those given in the sceneEye object. This is
%    necessary, since the optical image returned by piRender will have
%    parameters such as f-number, focal length, etc. set to default.
%    We want the values to reflect the retinal model we rendered.
%
% Inputs:
%    SE                    - The sceneEye object
%    ieObject              - The rendered optical image; typically
%                            from piRender.
%
% Outputs:
%   ieObject               - Object. The optical image is returned with the
%                            correct parameters set. 
%
% Optional key/value pairs:
%    meanIlluminancePerMM2 - The mean illuminance per mm^2 of pupil area
%                            Default 5 cd/m2. (Numeric)
%    scaleIlluminance      - Scale the returned illuminance with respect to
%                            the pupil area. Default true. 
%
% See also
%

%% Initialization

varargin = ieParamFormat(varargin);

p = inputParser;

p.addRequired('obj');
p.addRequired('ieObject');

p.addParameter('meanilluminancepermm2', 5, @isnumeric);
p.addParameter('scaleilluminance', true, @islogical);

p.parse(SE, ieObject, varargin{:});
meanIlluminancePerMM2 = p.Results.meanilluminancepermm2;
scaleIlluminance = p.Results.scaleilluminance;

%% Initial calculations
% ieObject = oiSet(ieObject, 'name', ...
%     sprintf('%s-%s', obj.name, datestr(now, 'mmm-dd, HH:MM')));

% Scene distance. Set it to infinity, since it doesn't apply to raytracing.
ieObject = oiSet(ieObject, 'distance', Inf);

% If there is a crop window we're going to have to adjust the FOV
crop_window   = SE.recipe.get('cropwindow');
retDistance   = SE.recipe.get('retina distance','m');
pupilDiameter = SE.recipe.get('pupil diameter','m');

full_size = 2 * tand(SE.recipe.get('fov') / 2) * retDistance;
image_size = full_size .* crop_window;

% Convert coordinates so center is (0, 0)
image_size = image_size-(full_size / 2);

% Assume image is square
field_angle(1) = atand(image_size(1) / retDistance);
field_angle(2) = atand(image_size(2) / retDistance);
fov_crop = abs(field_angle(1) - field_angle(2));

% For the optics focal length, we use the distance between the back of
% the lens and the retina. Although it is not exactly the same as the
% focal length for the eye (which also changes with accommodation), it
% is a good approximation for the oiWindow.
ieObject = oiSet(ieObject, 'optics focal length', retDistance);
ieObject = oiSet(ieObject, 'optics fnumber', retDistance / pupilDiameter);
ieObject = oiSet(ieObject, 'fov', fov_crop);

% Set optics that do not apply to the iset3d optical
% image. We may want to add these in in the future.
ieObject = oiSet(ieObject, 'optics model', 'iset3d');
ieObject = oiSet(ieObject, 'optics name', 'PBRT Navarro Eye');
ieObject = oiSet(ieObject, 'optics focal length',retDistance);

ieObject = oiSet(ieObject,'optics OTF', []);
ieObject = oiSet(ieObject,'optics name',SE.recipe.get('lens file'));
ieObject = oiSet(ieObject,'optics offaxis method','');
% ieObject = oiSet(ieObject,'optics vignetting',[]);

%% Calculate and apply lens transmittance

thisLens = Lens;
thisLens.wave    = oiGet(ieObject,'wave');
thisLens.density = SE.lensDensity;
ieObject         = oiSet(ieObject,'optics lens',thisLens);

%% Apply lens transmittance

irradiance    = oiGet(ieObject, 'photons');
transmittance = oiGet(ieObject, 'optics transmittance');

% We could conditionalize on any(transmittance(:) ~= 1)
transmittance = reshape(transmittance, [1 1 length(transmittance)]);
irradiance    = bsxfun(@times, irradiance, transmittance);
ieObject      = oiSet(ieObject, 'photons', irradiance);

%% Scale the irradiance with pupil size

% This is already done in piDat2ISET.m in some cases.
% 
% However for the human eye the lens file and other properties are not
% always available, and we are unable to extract the aperture (at the
% moment) in piDat2ISET.m. Until we we upgrade that code, we adjust the
% oi illuminance here because the pupil diameter is known in the sceneEye
% object. 
if(scaleIlluminance)
    lensArea = pi * (pupilDiameter*1e3 / 2) ^ 2;   % In mm^2
    meanIlluminance = meanIlluminancePerMM2 * lensArea;
    
    ieObject = oiAdjustIlluminance(ieObject, meanIlluminance);
    ieObject = oiSet(ieObject,'illuminance',oiCalculateIlluminance(ieObject));
end
% oiWindow(ieObject);

end
