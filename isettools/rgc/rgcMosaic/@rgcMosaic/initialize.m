function initialize(obj, rgc, cellTypeInd)
% intialize: a method of @rgcMosaic that initializes the object based on a
% series of input parameters that can include the location of the
% retinal patch.
% 
%             obj.initialize(rgc, scene, sensor, outersegment, varargin{:});
% 
% First, the name of the cell type is assigned based on the value passed in
% the type parameter. Next, spatial receptive fields are generated for the
% array of RGCs. Finally, the RGB temporal impulse responses for the center
% and surround are generated.
% 
% Inputs: 
%       scene: an isetbio scene structure
%       sensor: an isetbio sensor structure
%       os: an isetbio outer segment structure
%    Optional but recommended:
%       eyeSide: 'left' or 'right', which eye the retinal patch is from
%       patchRadius: radius of retinal patch in microns
%       patchAngle: polar angle of retinal patch
%     [These inputs determine the size of spatial receptive fields, and are
%       necessary to accurately model physiological responses.]
% 
% Outputs: the rgc object.
%  
% Example: from rgcGLM.m initiailize:
%        obj.mosaic{cellTypeInd} = rgcMosaicGLM(cellTypeInd, obj, scene, sensor, outersegment, varargin{:});
% 
% (c) isetbio
% 09/2015 JRG



namesCellTypes = {'onParasol';'offParasol';'onMidget';'offMidget';'smallBistratified'};
obj.cellType = namesCellTypes{cellTypeInd};

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

patchSizeX = rgc.spacing; %sensorGet(sensor, 'width', 'um');
sensorRows = rgc.row;% sensorGet(sensor,'rows');


umPerSensorPx = patchSizeX/sensorRows;

obj.rfDiameter = rfSizeMult(cellTypeInd)*(receptiveFieldDiameterParasol2STD/2)/umPerSensorPx; % in microns; divide by umPerScenePx to get pixels



% Build spatial RFs of all RGCs in this mosaic
[obj.sRFcenter, obj.sRFsurround, obj.rfDiaMagnitude, obj.cellLocation] = ...
    buildSpatialRFArray(rgc.spacing, rgc.row, rgc.col, obj.rfDiameter);


%% Add temporal response function to R, G, B channels
integrationTime = rgc.timing;%sensorGet(sensor, 'integration time');
for rgbInd = 1:3
    multFactor = rgbTempMult(rgbInd)*rfTempMult(cellTypeInd); % account for RGB channel and ON/OFF channel
    obj.tCenter{rgbInd,1} = multFactor*buildTemporalImpulseResponse(integrationTime);
    obj.tSurround{rgbInd,1} = multFactor*buildTemporalImpulseResponse(integrationTime);

end


