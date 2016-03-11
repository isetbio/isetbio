function initialize(obj, innerRetina, cellType)
%% Initialize an rgcMosaic for a particular cell type
%
%  DEPRECATED for rgcInitSpace/rgcInitTime
%
%       initialize(obj, innerRetina, cellType) 
%           [only called internally from rgcMosaic.m]
% 
% The object is intialized based on a series of input parameters that
% can include the location of the retinal patch.
% 
% First, the name of the cell type is assigned based on the value passed in
% the type parameter. Next, spatial receptive fields of the appropriate
% size are generated for the array of RGCs of that particular type. Then
% the RGB temporal impulse responses for the center and surround are
% generated.
% 
% Inputs: 
%       obj: the rgcMosaic object
%       innerRetina: the ir to which the rgcMosaic object is attached
%       cellType: determines the spatial RF and temporal impulse response
%           paramters of the rgcMosaic object; one of the following
%           strings:
% 
%              ON Parasol
%              OFF Parasol
%              ON Midget
%              OFF Midget
%              Small Bistratified
% 
% Outputs: the rgcMosaic object, which is then attached to the ir in
%       rgcMosaicCreate.
%  
% Example: [only called internally from rgcMosaic*]
%
%       innerRetina = rgcMosaicCreate(innerRetina,'mosaicType','on parasol');
%       innerRetina.mosaicCreate('mosaicType','on midget');
% 
% (c) isetbio
% 09/2015 JRG

%% Switch cell type string to index number 
% The index number helps with the generation of the receptive fields and
% impulse responses of the appropriate parameters for the cell type.
obj.cellType = cellType;
switch ieParamFormat(cellType)
    case{'onparasol'}
        cellTypeInd = 1;
    case{'offparasol'}        
        cellTypeInd = 2;
    case{'onmidget'}        
        cellTypeInd = 3;
    case{'offmidget'}
        cellTypeInd = 4;
    case{'smallbistratified'}
        cellTypeInd = 5;
    otherwise
        error('Unknown cell type');
end%switch


%% Generate spatial RFs of the approrpiate size for the cell type and TEE

% RF size scale parameters: Parasol RFs are the largest, while Midget RFs
% are about half the diameter of Parasol RFs, and SBC RFs are 1.2X the
% diameter of Parasol RFs.
% See Chichilnisky, E. J., and Rachel S. Kalmar. "Functional asymmetries 
% in ON and OFF ganglion cells of primate retina." The Journal of 
% Neuroscience 22.7 (2002).
% [ON Parasol; OFF Parasol; ON Midget; OFF Midget; Small bistratified];
rfSizeMult = [1 1 0.5 0.5 1.2];   % account for size differences between types
 
% Calcualte spatial RF diameter for ON Parasol cell at a particular TEE
% See Chichilnisky, E. J., and Rachel S. Kalmar. "Functional asymmetries 
% in ON and OFF ganglion cells of primate retina." The Journal of 
% Neuroscience 22.7 (2002), Fig. 5, pg. 2741.
receptiveFieldDiameterParasol2STD = receptiveFieldDiameterFromTEE(innerRetina.temporalEquivEcc);

% If os is osIdentity, determine number of pixels per spatial RF diameter.
% If os is osLinear or osBioPhys, determine the number of cones per spatial
% RF diameter.
patchSizeX = innerRetina.spacing; 
sensorRows = innerRetina.row;
umPerSensorPx = patchSizeX/sensorRows;

% The spatial RF diameter in pixels (or cones) is therefore the diameter in
% microns divided by the number of microns per pixel (or cone), scaled by
% the factor determined by the type of mosaic that is being created.
obj.rfDiameter = rfSizeMult(cellTypeInd)*(receptiveFieldDiameterParasol2STD/2); % in microns; divide by umPerScenePx to get pixels

% Build spatial RFs of all RGCs in this mosaic
[obj.sRFcenter, obj.sRFsurround, obj.rfDiaMagnitude, obj.cellLocation] = ...
    buildSpatialRFArray(innerRetina.spacing, innerRetina.row, innerRetina.col, obj.rfDiameter);


%% Generate temporal impulse response functions for R, G, B channels

% Temporal impulse response peak magnitude scale parameters.
% see "Spatial Properties and Functional Organization of Small
% Bistratified Ganglion Cells in Primate Retina", Field, et al.,
% J. Neuroscience, 2007, Fig. 1.

rfTempMult = [1 -1 1 -1 1];       % invert IR for OFF paraosl and midget

% The temopral impulse response function for each cell consists of three
% vectors, one for each RGB channel. They have similar shapes but different
% polarities and relative magnitudes depending on the cell type.
switch ieParamFormat(obj.cellType)
    case{'smallbistratified'}
        % SBCs have the B channel impulse response reversed in polarity
        rgbTempMult = [-0.4 -0.4 1]; % [R G B] magnitudes
    otherwise
        % The other four of the big five types have the same polarity
        rgbTempMult = [0.4 1 0.4]; % [R G B] magnitudes
end
       
% Get the sampling interval set in the outersegment or sensor
integrationTime = innerRetina.timing;
for rgbInd = 1:3
    % scale for differences in RGB channel and ON/OFF type
    multFactor = rgbTempMult(rgbInd)*rfTempMult(cellTypeInd); 
    % Build the separate impulse responses for center and surround; usually
    % the same.
    obj.tCenter{rgbInd,1} = multFactor*buildTemporalImpulseResponse(integrationTime);
    obj.tSurround{rgbInd,1} = multFactor*buildTemporalImpulseResponse(integrationTime);

end


