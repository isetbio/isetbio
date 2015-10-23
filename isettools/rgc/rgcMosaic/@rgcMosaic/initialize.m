function initialize(obj, type, rgc, scene, sensor, outersegment, varargin)
% intialize: a method of @rgcMosaic that initializes the object based on a
% series of input parameters that can include the location of the
% retinal patch.
% 
% First, the name of the cell type is assigned based on the value passed in
% the type parameter. Next, spatial receptive fields are generated for the
% array of RGCs. Finally, the RGB temporal impulse responses for the center
% and surround are generated.
% 
% Inputs: cell type index, rgc, sensor, outersegment.
% 
% Outputs: initialized rgc object.
% 
% Example:
% rgc1 = rgcLinear(sensor, osIdentity, 'right', 3.75, 180);
% rgc2 = rgcLNP(sensor, osIdentity, 'right', 3.75, 180);
% rgc3 = rgcGLM(sensor, osIdentity, 'right', 3.75, 180);
% 
% (c) isetbio
% 09/2015 JRG

%% Spatial RF arrays

namesCellTypes = {'onParasol';'offParasol';'onMidget';'offMidget';'smallBistratified'};
obj.cellType = namesCellTypes{type};

% Assign multipliers for size of spatial RFs and magnitudes of temporal IRs
% [ON Parasol; OFF Parasol; ON Midget; OFF Midget; Small bistratified];
rfSizeMult = [1 1 0.5 0.5 1.2];   % account for size differences between types
rfTempMult = [1 -1 1 -1 1];       % invert IR for OFF paraosl and midget
% rgbTempMult = [0.4 1 0.4];        % weight RGB components of temporal IR

% see "Spatial Properties and Functional Organization of Small
% Bistratified Ganglion Cells in Primate Retina", Field, et al.,
% J. Neuroscience, 2007, Fig. 1.
switch ieParamFormat(obj.cellType)
    case{'smallbistratified'}
        rgbTempMult = [-0.4 -0.4 1];
    otherwise
        rgbTempMult = [0.4 1 0.4];
end
        
% Calcualte spatial RF diameter
receptiveFieldDiameterParasol2STD = receptiveFieldDiameterFromTEE(rgc.temporalEquivEcc);

patchSizeX = sensorGet(sensor, 'width', 'um');
sceneRows = sceneGet(scene,'rows');
umPerScenePx = patchSizeX/sceneRows;

obj.rfDiameter = rfSizeMult(type)*(receptiveFieldDiameterParasol2STD/2)/umPerScenePx; % in microns; divide by umPerScenePx to get pixels



% Build spatial RFs of all RGCs in this mosaic
[obj.sRFcenter, obj.sRFsurround, obj.rfDiaMagnitude, obj.cellLocation] = ...
    buildSpatialRFArray(scene, sensor, obj.rfDiameter);


%% Add temporal response function to R, G, B channels
integrationTime = sensorGet(sensor, 'integration time');
for rgbInd = 1:3
    multFactor = rgbTempMult(rgbInd)*rfTempMult(type); % account for RGB channel and ON/OFF channel
    obj.tCenter{rgbInd,1} = multFactor*buildTemporalImpulseResponse(integrationTime);
    obj.tSurround{rgbInd,1} = multFactor*buildTemporalImpulseResponse(integrationTime);

end
