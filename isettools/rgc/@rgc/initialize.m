function initialize(obj, params)
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
% 09/2015 JRG Copyright ISETBIO Team

% p = inputParser;
% p.addRequired('outersegment');
% p.addParameter('eyeSide',  'left', @ischar);
% p.addParameter('eyeRadius', 1.25,  @isnumeric);
% p.addParameter('eyeAngle',  50,    @isnumeric);
% p.addParameter('name',     'macaque RGC', @ischar);

% p.parse(outersegment,varargin{:});
obj.eyeSide   = params.eyeSide;
obj.eyeRadius = params.eyeRadius;
obj.eyeAngle  = params.eyeAngle;
obj.name      = params.name;
obj.row       = params.row;
obj.col       = params.col;
obj.spacing   = params.spacing;
obj.timing    = params.timing;

% Give the object a name and slots for the five cell types
obj.mosaic = cell(1); % populated in initialize()

% Use the outersegment type to specify the inputs for the computation 
% obj.input = class(outersegment);
    
% Get the TEE.
obj.temporalEquivEcc = retinalLocationToTEE(obj.eyeAngle, obj.eyeRadius, obj.eyeSide);

%%

