function oi = opticsDLCompute(scene, oi)
%Diffraction limited optical image computation
%
%   oi = opticsDLCompute(scene,oi)
%
% The diffraction limited optical image calculation uses only a few
% parameters (f-number, focal length) to calculate the optical image.  The
% diffraction limited OTF is calculated on the fly in dlMTF, and applied to
% the scene image in this routine.
%
% See also:  oiCompute, opticsSICompute, opticsRayTrace
%
% Examples:
%  We don't normally call this function directly.  We use oiCompute which
%  then makes the appropriate call based on how the optics model is set.
%
% Copyright ImagEval Consultants, LLC, 2005

% TODO:  We should insert a geometric distortion function in this code,
% rather than using it only in the ray trace methods. 
if notDefined('scene'), scene = vcGetObject('scene'); end
if notDefined('oi'),    oi = vcGetObject('oi');       end
showWaitBar = ieSessionGet('waitbar');

opticsModel = oiGet(oi,'optics model');
if ~(strcmpi(opticsModel,'dlmtf') || strcmpi(opticsModel,'diffractionlimited'))
    error('Bad DL optics model %s',opticsModel);
end
wStr = 'OI-DL: '; 

% Set the basic parameters of the oi from the scene parameters.
oi = oiSet(oi,'wangular', sceneGet(scene,'wangular'));
oi = oiSet(oi,'data wave',sceneGet(scene,'wave'));

% Calculate the irradiance of the optical image in photons/(s m^2 nm)
if showWaitBar, wBar = waitbar(0,[wStr,' Calculating irradiance...']); end
oi   = oiSet(oi,'cphotons',oiCalculateIrradiance(scene,oi));

% We should insert a geometric distortion function here. This could be a
% button on the window that indicates we want the distortion computed.
% Which point in the computation should have the distortion? Before or
% after application of the OTF?
%
if showWaitBar
    waitbar(0.3,wBar,[wStr,' Calculating off-axis falloff']);
end

% Apply the offaxis fall-off.
%
% We either apply a standard cos4th calculation, or we use the more
% elaborate relative illumination derived from CodeV. stored inside of
% data\optics\Lens Design\Standard.
offaxismethod = oiGet(oi,'optics off axis method');
switch lower(offaxismethod)
    case {'skip','none',''}
    case 'cos4th'
        oi = opticsCos4th(oi);
    otherwise
        fprintf('Unknown offaxis method: %s. Using cos4th',optics.offaxis);
        oi = opticsCos4th(oi);
end

% We apply the MTF here.
if showWaitBar, waitbar(0.6,wBar,[wStr,' Applying OTF']); end
oi = opticsOTF(oi,scene);

% Diffuser and illuminance, or just illuminance.  Diffuser always resets
% the illuminance, which seems proper.
switch lower(oiGet(oi,'diffuserMethod'))
    case 'blur'
        if showWaitBar, waitbar(0.75,wBar,[wStr,' Diffuser']); end
        blur = oiGet(oi,'diffuserBlur','um');
        if ~isempty(blur), oi = oiDiffuser(oi,blur); end
    case 'birefringent'
        if showWaitBar, waitbar(0.75,wBar,[wStr,' Birefringent Diffuser']); end
        oi = oiBirefringentDiffuser(oi);
    case 'skip'
    otherwise
        error('Unknown diffuser method %s\n',oiGet(oi,'diffuser method'));
end

% Compute image illuminance (in lux)
% waitbar(0.9,wBar,[wStr,' Calculating illuminance']);
oi = oiSet(oi,'illuminance',oiCalculateIlluminance(oi));

if showWaitBar, close(wBar); end

end