function initialize(obj, type, rgc, sensor, outersegment, varargin)
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

%% Spatial RF arrays

obj.nameCellType = rgc.namesCellTypes{type};

% [ON Parasol; OFF Parasol; ON Midget; OFF Midget; Small bistratified];
rfSizeMult = [1 1 0.5 0.5 1.2];   % account for size differences between types
rfTempMult = [1 -1 1 -1 1];       % invert IR for OFF paraosl and midget
rgbTempMult = [0.4 1 0.4];        % weight RGB components of temporal IR
        
receptiveFieldDiameterParasol2STD = receptiveFieldDiameterFromTEE(rgc.temporalEquivEcc);
obj.receptiveFieldDiameter1STD = rfSizeMult(type)*receptiveFieldDiameterParasol2STD/2;

[obj.spatialRFArray obj.spatialRFonedim ...
    obj.spatialRFcontours obj.spatialRFFill obj.cellCenterLocations] = ...
    buildSpatialRFArray(sensor, obj.receptiveFieldDiameter1STD);


%% Add temporal response function to R, G, B channels
integrationTime = sensorGet(sensor, 'integration time');
for rgbInd = 1:3
    multFactor = rgbTempMult(rgbInd)*rfTempMult(type); % account for RGB channel and ON/OFF channel
    obj.temporalImpulseResponseCenterRGB{rgbInd,1} = multFactor*buildTemporalImpulseResponse(integrationTime);
    obj.temporalImpulseResponseSurroundRGB{rgbInd,1} = multFactor*buildTemporalImpulseResponse(integrationTime);

end
