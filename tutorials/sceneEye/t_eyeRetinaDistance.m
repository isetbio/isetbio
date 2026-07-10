%% t_eyeRetinaDistance.m
% SkipFile
% Human-eye retinal-distance sweep currently relies on legacy sceneEye/
% docker setup and many CPU optics renders, so keep it out of smoke tests.
%
% This tutorial renders a retinal image of "slanted edge." We can use
% this slanted bar to estimate the modulation transfer function of the
% optical system.  But in this script we mainly explore what happens
% as we set the accommodation and retinal distance of the human eye
% models (Navarro and Arizona).
%
% Also of note, the fringe color changes nicely as we sweep out
% different retinal distances.
%
% It is worth noting that small changes in the retinal distance (50
% microns) have substantial effects on the fringing. The idea that the
% visual system can absolutely count on the fringing for a precise
% measurement, or that different people are the same, seems unlikely
% to me (BW).
%
% Depends on: ISETBIO, Docker, ISETCam
%
%  
% See also
%   t_eyeArizona, t_eyeNavarro
%

%% Check ISETBIO and initialize

ieInit;
if ~piDockerExists, piDockerConfig; end

%% Set up the slanted bar scene

% This can be set for all of the eye models.
modelName = {'navarro','arizona','legrand'};
mm = 3;

thisSE = sceneEye('slantedEdge','eye model',modelName{mm});
thisSE.set('to',[0 0 0]);

thisLight = piLightCreate('spot light 1', 'type','spot','rgb spd',[1 1 1]);
thisSE.set('light',thisLight, 'add');
thisSE.set('light',thisLight.name,'specscale',0.5);

% Set up the image
thisSE.set('fov',1);                % Field of view
thisSE.set('spatial samples',[256, 256]);  % Number of OI sample points
thisSE.set('film diagonal',2);          % mm
thisSE.set('rays per pixel',256);
thisSE.set('n bounces',2);

thisSE.set('lens density',0);       % Remove pigment. Yellow irradiance is harder to see.
thisSE.set('diffraction',false);
thisSE.set('pupil diameter',3);

% We run this for a distance edge.  We set the retina distance - to
% have the right aberration for a far away edge.
oDistance = 10;
thisSE.set('object distance',oDistance);
fprintf('Object distance %.2f m\n',oDistance);

%{
 piAssetGeometry(thisSE.recipe);
 thisSE.recipe.show('lights');
%}

%% Render with model eye, varying diffraction setting

% Use model eye
thisSE.set('use optics',true);

thisSE.set('fov',1);                % Small Field of view

% This sets the chromaticAberrationEnabled flag.
% Works in V4 - May 28, 2023 (ZL)
nSpectralBands = 8;
thisSE.set('chromatic aberration',nSpectralBands);

if mm < 3, thisSE.set('accommodation',1/oDistance); end
thisSE.get('focal distance')

% thisSE.set('retina distance',16.32);  % Default
%{
Navarro eye Obj 10m the in focus retina distance is 16.35 - 16.4mm.
Arizona eye Obj 10m the in focus retina distance is 16.55mm.
Legrand eye Obj 10m the in focus retina distance is 16.35 mm
%}

% We sweep the retinal distance by 50 micron steps to see the blur
% change.
thisDocker = isetdocker;
for rr =  16.1:0.05:16.75
    thisSE.set('retina distance',rr);
    name = sprintf('%s Foc %.2f Obj %.2f Ret %0.2f',modelName{mm}(1:2),...
        thisSE.get('focal distance'),...
        oDistance,...
        thisSE.get('retina distance','mm'));
    thisSE.summary;
    oi = thisSE.piWRS('name',name,'docker',thisDocker,'show',true);
end

% If you want to reduce the rendering noise:
%
%   oi = ieGetObject('oi'); oi = piAIdenoise(oi); ,oiWindow(oi);
%

%% END
