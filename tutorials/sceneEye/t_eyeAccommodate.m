%% t_eyeAccommodate.m
% SkipFile
% Uses obsolete dockerWrapper human-eye setup; needs sceneEye rewrite.
%
% This tutorial renders a retinal image of "slanted edge" to check how
% close the accommodation parameters match
%
% We use the slanted bar set at 1 m from the eye.  We adjust the eye
% model accommodation to bring the 1 m slanted edge into good focus.
%
% The setting for the 'best' focus appears to be not 1 Diopter (e.g.,
% 1m away) but closer to 1.2 D.
%
% The Navarro model already has a conversion factor for this purpose.
% The Arizona model does not.  I am considering whether or not to
% implement an updated conversion for both of these models (BW).
%
% The LeGrand eye model does not have an accommodation.  It only works
% for objects at a distance.
%
% Depends on: ISETBIO, Docker, ISETCam
%
% See also
%   t_eyeRetinaDistance, t_eyeArizona, t_eyeNavarro
%

%% Check ISETBIO and initialize

ieInit;
if ~piDockerExists, piDockerConfig; end

%% Set up scene and eye model

% Only these eye models (not legrand) can accommodate
modelName = {'navarro','arizona'};
mm = 1;

% Choose 1 or 2 for Navarro or Arizona
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

% We run this for a closer distance to see how close the accommodation
% is to solving the new focal distance (object distance in focus).
oDistance = 1;
thisSE.set('object distance',oDistance);
fprintf('Object distance %.2f m\n',oDistance);

%{
thisSE.piWRS;
%}
%{
 piAssetGeometry(thisSE.recipe);
 thisSE.recipe.show('lights');
%}

%% Render with model eye, varying diffraction setting

% Use model eye
thisSE.set('use optics',true);

thisSE.set('fov',1);                % Field of view

% This sets the chromaticAberrationEnabled flag.
% Works in V4 - May 28, 2023 (ZL)
nSpectralBands = 8;
thisSE.set('chromatic aberration',nSpectralBands);

thisSE.set('accommodation',1/oDistance);
thisSE.get('focal distance')
inFocusAcc = 1/oDistance;   % 1 diopter

%{
 Navarro eye, Obj at 1 diopter, best accommodation is 1.2 D
 Arizona eye, Obj at 1 diopter, best accommodation is 1.25 D.
%}

% We step the accommodation, hoping that the best focus matches the
% object distance (1 diopter)
thisDocker = isetdocker;
for aa =  0.8:.2:1.2
    thisSE.set('accommodation',aa);
    name = sprintf('%s Foc %.2f Obj %.2f (D)',modelName{mm}(1:2),...
        thisSE.get('accommodation'));
    thisSE.summary;
    oi = thisSE.piWRS('name',name,'docker',thisDocker,'show',true);
end

% oi = ieGetObject('oi'); oi = piAIdenoise(oi); ,oiWindow(oi);


%% END
