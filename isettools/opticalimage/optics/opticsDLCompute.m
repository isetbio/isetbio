function oi = opticsDLCompute(scene, oi)
% Diffraction limited optical image computation
%
% Syntax:
%   oi = opticsDLCompute(scene, oi)
%
% Description:
%    The diffraction limited optical image calculation uses only a few
%    parameters (f-number, focal length) to calculate the optical image.
%    The diffraction limited OTF is calculated on the fly in dlMTF, and
%    applied to the scene image in this routine.
%
% Inputs:
%    scene - Struct. A scene structure
%    oi    - Struct. An optical image structure
%
% Outputs:
%    oi    - Struct. The modified optical image structure
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * We don't normally call this function directly. We use oiCompute
%      which then makes the appropriate call based on how the optics model
%      is set.
%    * TODO: Insert geometric distortion function rather than restricting
%      it to the ray trace methods.
%
% See Also:
%    oiCompute, opticsSICompute
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005
%    03/12/18  jnm  Formatting


if notDefined('scene'), scene = vcGetObject('scene'); end
if notDefined('oi'), oi = vcGetObject('oi'); end
showWaitBar = ieSessionGet('waitbar');

opticsModel = oiGet(oi, 'optics model');
if ~(strcmpi(opticsModel, 'dlmtf') || ...
        strcmpi(opticsModel, 'diffractionlimited'))
    error('Bad DL optics model %s', opticsModel);
end
wStr = 'OI-DL: ';

% Set the basic parameters of the oi from the scene parameters.
oi = oiSet(oi, 'wangular', sceneGet(scene, 'wangular'));
oi = oiSet(oi, 'data wave', sceneGet(scene, 'wave'));

% Calculate the irradiance of the optical image in photons/(s m^2 nm)
if showWaitBar
    wBar = waitbar(0, [wStr, ' Calculating irradiance...']);
end
oi   = oiSet(oi, 'photons', oiCalculateIrradiance(scene, oi));

% [Note: XXX - We should insert a geometric distortion function here. This
% could be a button on the window that indicates we want the distortion
% computed. Which point in the computation should have the distortion?
% Before or after application of the OTF?]
if showWaitBar
    waitbar(0.3, wBar, [wStr, ' Calculating off-axis falloff']);
end

% Apply the offaxis fall-off.
%
% We either apply a standard cos4th calculation, or we use the more
% elaborate relative illumination derived from CodeV. stored inside of
% data\optics\Lens Design\Standard.
offaxismethod = oiGet(oi, 'optics off axis method');
switch lower(offaxismethod)
    case {'skip', 'none', ''}
    case 'cos4th'
        oi = opticsCos4th(oi);
    otherwise
        fprintf('Unknown offaxis method: %s. Using cos4th', ...
            optics.offaxis);
        oi = opticsCos4th(oi);
end

% We apply the MTF here.
if showWaitBar, waitbar(0.6, wBar, [wStr, ' Applying OTF']); end
oi = opticsOTF(oi, scene);

% Diffuser and illuminance, or just illuminance. Diffuser always resets the
% illuminance, which seems proper.
switch lower(oiGet(oi, 'diffuserMethod'))
    case 'blur'
        if showWaitBar, waitbar(0.75, wBar, [wStr, ' Diffuser']); end
        blur = oiGet(oi, 'diffuserBlur', 'um');
        if ~isempty(blur), oi = oiDiffuser(oi, blur); end
    case 'birefringent'
        if showWaitBar
            waitbar(0.75, wBar, [wStr, ' Birefringent Diffuser']);
        end
        oi = oiBirefringentDiffuser(oi);
    case 'skip'
    otherwise
        error('Unknown diffuser method %s\n', ...
            oiGet(oi, 'diffuser method'));
end

% Compute image illuminance (in lux)
% waitbar(0.9, wBar, [wStr, ' Calculating illuminance']);
oi = oiSet(oi, 'illuminance', oiCalculateIlluminance(oi));

if showWaitBar, close(wBar); end

end