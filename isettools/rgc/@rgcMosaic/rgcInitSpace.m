function rgcM = rgcInitSpace(rgcM,rgcLayer,cellType,varargin)
% Initialize the spatial rf properties of a rgc mosaic for a cell type
%
%    rgcM = rgcInitSpace(rgcM,rgcLayer,cellType)
%
% cellType is one of (spacing and case are irrelevant)
%  {'ON Parasol', 'OFF Parasol', 'ON Midget', OFF Midget', 'Small bistratified'}
%
% RF size scale parameters: Parasol RFs are the largest, while Midget RFs
% are about half the diameter of Parasol RFs, and SBC RFs are 1.2X the
% diameter of Parasol RFs.
%
% In order to create the spatial RFs for each type of cell, we compute the
% RF size for the On Parasol cells at a particular eccentricity based on
% data from Watanabe & Rodieck (J. Comp. Neurol., 1989), Croner & Kaplan
% (Vision Research, 1995) and Chichilniksy & Kalmar (J Neurosci., 2002,
% Fig. 5, pg. 2741.).
%
% For other types of RGCs, we multiply the On Parasol RF size by a unitless
% factor as listed here. The multipliers are based on Dacey, 2004, "Retinal
% ganglion cell diversity", in The Cognitive Neurosciences, as well as
% parameter fits to mosaic data from the Chichilnisky Lab found in the
% EJLExperimentalRGC iestbio repository.
%
% See Chichilnisky, E. J., and Rachel S. Kalmar. "Functional asymmetries
% in ON and OFF ganglion cells of primate retina." The Journal of
% Neuroscience 22.7 (2002).
%
% Example:
% See also:
%
% JRG/BW ISETBIO Team 2015

%% All three arguments are required

p = inputParser;
p.KeepUnmatched = true;
vFunc = @(x)(ismember(class(x),{'rgcGLM','rgcLNP'}));
p.addRequired('rgcM',vFunc);
p.addRequired('rgcLayer',@(x)(isequal(class(x),'rgcLayer')));

vFunc = @(x)(ismember(ieParamFormat(x),...
    {'onparasol', 'offparasol', 'onmidget', 'offmidget', 'smallbistratified'}));
p.addRequired('cellType',vFunc);
p.addParameter('inMosaic',1,@isscalar);
p.addParameter('rfDiameter',[],@isnumeric); % microns
p.parse(rgcM,rgcLayer,cellType,varargin{:});

rgcM.rfDiameter = p.Results.rfDiameter;
inMosaic = p.Results.inMosaic;

%% Set up defaults for the sizes and weights.
switch ieParamFormat(cellType)
    case 'onparasol'
        rfSizeMult = 1;      % On Parasol RF size multiplier
    case 'offparasol'
        rfSizeMult = 0.85;   % Off Parasol RF size multiplier
    case 'onmidget'
        rfSizeMult = 0.5;    % On Midget RF size multiplier
    case 'offmidget'
        rfSizeMult = 0.425;  % Off Midget RF size multiplier
    case 'smallbistratified'
        rfSizeMult = 1.1;    % SBC RF size multiplier
end

% If diameter has not been set, assign its value based on the literature
if isempty(rgcM.rfDiameter)
    % Calculate spatial RF diameter for ON Parasol cell at a particular TEE
    % See Chichilnisky, E. J., and Rachel S. Kalmar. "Functional asymmetries
    % in ON and OFF ganglion cells of primate retina." The Journal of
    % Neuroscience 22.7 (2002), Fig. 5, pg. 2741. 2STD fit in micrometers.
    retinalLocationToTEE = temporalEquivEcc(rgcLayer.center);
    receptiveFieldDiameterParasol2STD = ...
        receptiveFieldDiameterFromTEE(retinalLocationToTEE);
    
    % The spatial RF diameter in pixels (or cones) is therefore the diameter in
    % microns divided by the number of microns per pixel (or cone), scaled by
    % the factor determined by the type of mosaic that is being created.
    
    % in micrometers; divide by umPerScenePx to get pixels
    rgcM.rfDiameter = rfSizeMult*(receptiveFieldDiameterParasol2STD/2);
end

% Build spatial RFs of all RGCs in this mosaic
% We need to know the number of rows and columns in the bipolarMosaic that
% drives this mosaic.  So, we need to know which bipolar mosaic we are
% using as input here.
%
% RIGHT NOW I JUST USE THE FIRST MOSAIC.  BUT THIS SHOULD GET PASSED IN AND
% SET!!!
sz = size(rgcLayer.input.mosaic{inMosaic}.cellLocation);
[rgcM.sRFcenter, rgcM.sRFsurround, rgcM.cellLocation, rgcM.tonicDrive, rgcM.ellipseMatrix] = ...
    buildSpatialRFArray(rgcLayer.size(2), ...
    sz(1), ...
    sz(2), ...
    rgcM.rfDiameter, ...
    p.Unmatched);

% sRFcenter, sRFsurround and cellLocation are in units of the inputObj - scene pixels,
% cones or bipolar cells.
% rfDiaMagnitude is in units of micrometers.
% tonicDrive is in units of conditional intensity.

end
