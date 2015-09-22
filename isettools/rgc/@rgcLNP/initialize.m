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
        
for cellTypeInd = 1:obj.numberCellTypes
    obj.mosaic{cellTypeInd,1} = struct(...
        'nameCellType',obj.namesCellTypes{cellTypeInd},...
        'receptiveFieldDiameter1STD',[],...
        'spatialRFArray',[],...
        'spatialRFonedim',[],...
        'spatialRFcontours',[],...
        'cellCenterLocations',[],...
        'temporalImpulseResponse',[],...
        'generatorFunction',[],...
        'linearResponse',[],...
        'nlResponse',[],...
        'spikeResponse',[]...
        );
end

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


%% Find the appropriate RF size, the 2 STD diameter in um.
receptiveFieldDiameterParasol2STD = receptiveFieldDiameterFromTEE(temporalEquivEcc);

receptiveFieldDiameterParasol1STD = receptiveFieldDiameterParasol2STD/2;

obj.mosaic{1}.receptiveFieldDiameter1STD = receptiveFieldDiameterParasol1STD;
obj.mosaic{2}.receptiveFieldDiameter1STD = receptiveFieldDiameterParasol1STD;

%% Hack for getting spatial RF arrays of parasol RGCs
% Still need to make ON/OFF magnitude reversal in builSpatialRFArray
rfMosaic = figure;
[obj.mosaic{1}.spatialRFArray obj.mosaic{1}.spatialRFonedim obj.mosaic{1}.spatialRFcontours obj.mosaic{1}.spatialRFFill obj.mosaic{1}.cellCenterLocations] = buildSpatialRFArray(sensor, receptiveFieldDiameterParasol1STD);
[obj.mosaic{2}.spatialRFArray obj.mosaic{2}.spatialRFonedim obj.mosaic{2}.spatialRFcontours obj.mosaic{2}.spatialRFFill obj.mosaic{2}.cellCenterLocations] = buildSpatialRFArray(sensor, receptiveFieldDiameterParasol1STD);

%% Hack for getting spatial RF arrays of midget RGCs
receptiveFieldDiameterMidget1STD = receptiveFieldDiameterParasol1STD/2;

obj.mosaic{3}.receptiveFieldDiameter1STD = receptiveFieldDiameterMidget1STD;
obj.mosaic{4}.receptiveFieldDiameter1STD = receptiveFieldDiameterMidget1STD;

[obj.mosaic{3}.spatialRFArray obj.mosaic{3}.spatialRFonedim obj.mosaic{3}.spatialRFcontours obj.mosaic{3}.spatialRFFill obj.mosaic{3}.cellCenterLocations] = buildSpatialRFArray(sensor, receptiveFieldDiameterMidget1STD);
[obj.mosaic{4}.spatialRFArray obj.mosaic{4}.spatialRFonedim obj.mosaic{4}.spatialRFcontours obj.mosaic{4}.spatialRFFill obj.mosaic{4}.cellCenterLocations] = buildSpatialRFArray(sensor, receptiveFieldDiameterMidget1STD);

%% Spatial RF for small bistratified

receptiveFieldDiameterSBC1STD = 2*receptiveFieldDiameterParasol1STD;

obj.mosaic{5}.receptiveFieldDiameter1STD = receptiveFieldDiameterSBC1STD;
[obj.mosaic{5}.spatialRFArray obj.mosaic{5}.spatialRFonedim obj.mosaic{5}.spatialRFcontours obj.mosaic{5}.spatialRFFill obj.mosaic{5}.cellCenterLocations] = buildSpatialRFArray(sensor, receptiveFieldDiameterSBC1STD);

%% Should loop through types and call function
% Still need to adjust input parameters to function

% for cellTypeInd = 1:obj.numberCellTypes
%     [obj.mosaic{cellTypeInd}.spatialRFArray obj.mosaic{cellTypeInd,1}.cellCenterLocations]]= buildSpatialRFArray(sensor, obj.receptiveFieldDiameter1STD, obj.namesCellTypes{cellTypeInd});
% end

%% Add temporal response function to R, G, B channels
tempSign = [1,-1,1,-1,1]; % invert response for OFF paraosl and midget
for cellTypeInd = 1:obj.numberCellTypes
    obj.mosaic{cellTypeInd,1}.temporalImpulseResponse = tempSign(cellTypeInd)*buildTemporalImpulseResponse(sensor, obj.namesCellTypes{cellTypeInd});
end

%% Add generator function
% Need to make this into Gaussian CDF
for cellTypeInd = 1:obj.numberCellTypes
    obj.mosaic{cellTypeInd,1}.generatorFunction = @exp;
end


