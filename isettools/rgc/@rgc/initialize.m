function initialize(obj, outersegment, varargin)
% iinitializes the rgc object. 
%
% Inputs: 
%   outersegment - outer segment object
%   eye parameters are sent in as key/value pairs
%     name, eyeSide , eyeRadius (in mm), eyeAngle (degrees)
%
% Outputs: initialized rgc object.
% 
% The user inputs the location of the retinal patch with (eye side, reitnal
% patch radius, patch angle), and the temporal equivalent eccentricity
% (TEE) is calculated from Chichilnisky & Kalmar, 2002, J. Neurosci, where
% 
%   TEE = sqrt((0.61*X^2)+Y^2) (corrected from the text of the paper)
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
% center and surround temporal impulse responses, and for LNP objects, the
% generator function, and for GLM objects, the post-spike and coupling
% filters.
% 
% 
% Example:
%  rgc = rgcCreate() ....
%  rgc = rgcLinear(osIdentity, 'eyeSide','right', 'eyeRadius',3.75, 'eyeAngle',180);
%  rgc = rgcLNP(sensor, osIdentity, 'right', 3.75, 180);
%  rgc = rgcGLM(sensor, osIdentity, 'right', 3.75, 180);
% 
% 09/2015 JRG


p = inputParser;
p.addRequired('outersegment');
p.addParameter('eyeSide',  'left', @ischar);
p.addParameter('eyeRadius', 1.25,  @isnumeric);
p.addParameter('eyeAngle',  50,    @isnumeric);
p.addParameter('name',     'macaque RGC', @ischar);

p.parse(outersegment,varargin{:});
eyeSide   = p.Results.eyeSide;
eyeRadius = p.Results.eyeRadius;
eyeAngle  = p.Results.eyeAngle;
obj.name  = p.Results.name;

% Give the object a name and slots for the five cell types
obj.mosaic = cell(1); % populated in initialize()

% Use the outersegment type to specify the inputs for the computation 
if isa(outersegment,'osIdentity'),  obj.input = 'rgb';
else                                obj.input = 'cone current';
end

% coneSize = sensorGet(sensor, 'pixel size', 'um' );
% patchSizeX = sensorGet(sensor, 'width', 'um');
% patchSizeY = sensorGet(sensor, 'height', 'um');
% fov = sensorGet(sensor,'fov');
% numCones = sensorGet(sensor, 'size');

%% Parasol ON or OFF RGCs.

% % Specify location of the retinal patch.
% if nargin < 3 % no user input, set up defaults
%     eyeSide = 'left';
%     eyeRadius  = 1.25; % in um
%     eyeAngle   = 50; % in degrees
% else    
%     eyeSide = varargin{1,1};
%     eyeRadius  = varargin{1,2}; % in um
%     eyeAngle   = varargin{1,3}; % in degrees
% end

% Set these as properties of the object.
% obj.side = leftOrRightEye;
    
% Get the TEE.
obj.temporalEquivEcc = retinalLocationToTEE(eyeAngle, eyeRadius, eyeSide);
    
% Plot the TEE and the location of the retinal patch.
% plotPatchEccentricity(retinalTheta, retinalRadius, leftOrRightEye, temporalEquivEcc)

%%

