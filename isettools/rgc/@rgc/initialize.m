function initialize(obj, scene, sensor, outersegment, varargin)
% initialize: a method of @rgc that initializes the object. The user inputs 
% the location of the retinal patch with (eye side, reitnal patch radius, 
% patch angle), and the temporal equivalent eccentricity (TEE) is calculated
% from Chichilnisky & Kalmar, 2002, J. Neurosci, where 
% 
% TEE = sqrt((0.61*X^2)+Y^2) (corrected from the text of the paper)
% 
% The TEE is used to cacluate the spatial receptive field (RF) diameter
% from Fig. 5 of Chichilnisky & Kalmar, 2002. This intiailization procedure
% calls the rgcMosaic initialization procedure to build the mosaics for
% five cell types:
% 
% 1. ON parasol
% 2. OFF parasol
% 3. ON midget
% 4. OFF midget
% 5. small bistratified
% 
% The mosaic object is initialized with the center and surround RF, the
% center and surround temporal impulse responses, and for LNP
% objects, the generator function, and for GLM objects, the post-spike and
% coupling filters.
% 
% Inputs: sensor, outersegment, eye side, retinal patch radius, patch angle.
% 
% Outputs: initialized rgc object.
% 
% Example:
% rgc1 = rgcLinear(sensor, osIdentity, 'right', 3.75, 180);
% rgc2 = rgcLNP(sensor, osIdentity, 'right', 3.75, 180);
% rgc3 = rgcGLM(sensor, osIdentity, 'right', 3.75, 180);
% 
% 09/2015 JRG

if isa(outersegment,'osIdentity')
    obj.input = 'rgb';
else 
    obj.input = 'cone current';
end

obj.name = 'macaque RGC';
obj.mosaic = cell(5,1); % populated in initialize()

% coneSize = sensorGet(sensor, 'pixel size', 'um' );
% patchSizeX = sensorGet(sensor, 'width', 'um');
% patchSizeY = sensorGet(sensor, 'height', 'um');
% fov = sensorGet(sensor,'fov');
% numCones = sensorGet(sensor, 'size');

%% Parasol ON or OFF RGCs.

% Specify location of the retinal patch.
if nargin < 3 % no user input, set up defaults
    leftOrRightEye = 'left';
    retinalRadius = 1.25; % in um
    retinalTheta = 50; % in degrees
else    
    leftOrRightEye = varargin{1,1};
    retinalRadius = varargin{1,2}; % in um
    retinalTheta = varargin{1,3}; % in degrees
end

% Set these as properties of the object.
% obj.side = leftOrRightEye;
    
% Get the TEE.
temporalEquivEcc = retinalLocationToTEE(retinalTheta, retinalRadius, leftOrRightEye);

obj.temporalEquivEcc = temporalEquivEcc;
    
% Plot the TEE and the location of the retinal patch.
plotPatchEccentricity(retinalTheta, retinalRadius, leftOrRightEye, temporalEquivEcc)

%%

