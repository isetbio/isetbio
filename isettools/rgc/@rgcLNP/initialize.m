function initialize(obj, sensor, outersegment, varargin)
% intialize: a method of @rgcLNP that initializes the object based on a
% series of input parameters that can include the location of the
% retinal patch.
% 
% Inputs:
% 
% Outputs:
% 
% Example:
% 
% (c) isetbio
% 09/2015 JRG

% mosaic = struct('onParasol',[],'offParasol',[],'onMidget',[],'offMidget',[],'smallBistratified',[]);
obj.animal = 'macaque';
obj.numberCellTypes = 5;
obj.namesCellTypes = {'onParasol';'offParasol';'onMidget';'offMidget';'smallBistratified'};
obj.mosaic = cell(obj.numberCellTypes,1); % populated in initialize()
        
coneSize = sensorGet(sensor, 'pixel size', 'um' );
patchSizeX = sensorGet(sensor, 'width', 'um');
patchSizeY = sensorGet(sensor, 'height', 'um');
fov = sensorGet(sensor,'fov');
numCones = sensorGet(sensor, 'size');

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
obj.eyeLeftOrRight = leftOrRightEye;
obj.patchLocationPolarRadiusMicrometers = retinalRadius;
obj.patchLocationPolarAngleDegrees = retinalTheta;
    
% Get the TEE.
temporalEquivEcc = retinalLocationToTEE(retinalTheta, retinalRadius, leftOrRightEye);

obj.temporalEquivEcc = temporalEquivEcc;
    
% Plot the TEE and the location of the retinal patch.
plotPatchEccentricity(retinalTheta, retinalRadius, leftOrRightEye, temporalEquivEcc)

%%



for cellTypeInd = 1:obj.numberCellTypes
    obj.mosaic{cellTypeInd} = rgcMosaicLNP(cellTypeInd, obj, sensor, outersegment, varargin{:});
end


