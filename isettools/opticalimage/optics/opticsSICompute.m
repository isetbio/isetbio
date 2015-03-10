function oi = opticsSICompute(scene,oi)
%Calculate OI irradiance using a custom shift-invariant PSF
%
%    oi = opticsSICompute(scene,oi)
%
% See also: opticsRayTrace, oiCompute, opticsOTF
%
% Example
%    scene = vcGetObject('scene');
%    oi    = vcGetObject('oi');
%    
% Copyright ImagEval Consultants, LLC, 2005

if notDefined('scene'), error('Scene required.'); end
if notDefined('oi'), error('Opticalimage required.'); end
showWbar = ieSessionGet('waitbar');

% This is the default compute path
optics = oiGet(oi,'optics');

% Compute the basic parameters of the oi from the scene parameters.
oi = oiSet(oi,'wangular',sceneGet(scene,'wangular'));
oi = oiSet(oi,'spectrum',sceneGet(scene,'spectrum'));

if isempty(opticsGet(optics,'otfdata')), error('No psf data'); end

% We use the custom data.
optics = opticsSet(optics,'spectrum',oiGet(oi,'spectrum'));
oi     = oiSet(oi,'optics',optics);

% Convert radiance units to optical image irradiance (photons/(s m^2 nm))
if showWbar
    wBar = waitbar(0,'OI-SI: Calculating irradiance...');
end

oi = oiSet(oi,'cphotons',oiCalculateIrradiance(scene,optics));

%-------------------------------
% Distortion would go here. If we included it.
%-------------------------------

if showWbar, waitbar(0.3,wBar,'OI-SI Calculating off-axis falloff'); end

% Now apply the relative illumination (offaxis) fall-off
% We either apply a standard cos4th calculation, or we skip.
%waitbar(0.3,wBar,'OI-SI: Calculating off-axis falloff');
offaxismethod = opticsGet(optics,'offaxismethod');
switch lower(offaxismethod)
    case {'skip','none',''}
    case 'cos4th'
        oi = opticsCos4th(oi);
    otherwise
        fprintf('Unknown offaxis method: %s\nUsing cos4th',optics.offaxis);
        oi = opticsCos4th(oi);
end

if showWbar, waitbar(0.6,wBar,'OI-SI: Applying OTF'); end
% This section applys the OTF to the scene radiance data to create the
% irradiance data.
%
% If there is a depth plane in the scene, we also blur that and put the
% 'blurred' depth plane in the oi structure.
if showWbar, waitbar(0.6,wBar,'Applying OTF-SI'); end
oi = opticsOTF(oi,scene);

switch lower(oiGet(oi,'diffuserMethod'))
    case 'blur'
       if showWbar, waitbar(0.75,wBar,'Diffuser'); end
        blur = oiGet(oi,'diffuserBlur','um');
        if ~isempty(blur), oi = oiDiffuser(oi,blur); end
    case 'birefringent'
       if showWbar, waitbar(0.75,wBar,'Birefringent Diffuser'); end
        oi = oiBirefringentDiffuser(oi);
    case 'skip'
        
end

% Compute image illuminance (in lux)
if showWbar, waitbar(0.9,wBar,'OI: Calculating illuminance'); end
oi = oiSet(oi,'illuminance',oiCalculateIlluminance(oi));

if showWbar, delete(wBar); end

end