function initialize(obj, os, params)
%% Initializes the ir object. 
%
% Inputs: 
%   outersegment - outer segment object
%   eye parameters are sent in as key/value pairs
%     name, eyeSide , eyeRadius (in mm), eyeAngle (degrees)
%
% Outputs: initialized ir object.
% 
% The user inputs the location of the retinal patch with (eye side, reitnal
% patch radius, patch angle), and the temporal equivalent eccentricity
% (TEE) is calculated from Chichilnisky & Kalmar, 2002, J. Neurosci, where
% 
%   TEE = sqrt((0.61*X^2)+Y^2) (corrected from the text of the paper)
% 
% The TEE is used to cacluate the spatial receptive field (RF) diameter
% from Fig. 5 of Chichilnisky & Kalmar, 2002. 
% 
% Mosaics are not created with the initialization of the ir object. They
% are created manually with the mosaicCreate method. See t_rgc.m.
% 
% Example:
% 
% 
%   os  = osCreate('identity');
%   innerRetina1 = irCreate(os,'GLM','name','myRGC'); 
% 
%   params.name = 'Macaque inner retina 1'; % 
%   innerRetina2 = ir(os, params);
% 
% 09/2015 JRG Copyright ISETBIO Team

%% Set properties
% Set properties from the input parser from irCreate
obj.eyeSide   = params.eyeSide;
obj.eyeRadius = params.eyeRadius;
obj.eyeAngle  = params.eyeAngle;
obj.name      = params.name;

% Set properties dependent on os object.
obj.spacing = osGet(os,'patch size'); % Cone width
obj.timing  = osGet(os,'time step'); % Temporal sampling
    
switch class(os)
    case{'osIdentity'}
        [obj.row, obj.col, ~, ~] = size(osGet(os,'rgbData'));
    otherwise    
        [obj.row, obj.col, ~] = size(osGet(os,'coneCurrentSignal'));
end

% Initialize the mosaic property but do not generate any mosaics
obj.mosaic = cell(1); % populated by mosaicCreate method
    
% Get the TEE.
obj.temporalEquivEcc = retinalLocationToTEE(obj.eyeAngle, obj.eyeRadius, obj.eyeSide);

%%

