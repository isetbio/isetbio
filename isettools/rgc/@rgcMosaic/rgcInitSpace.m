function rgcM = rgcInitSpace(rgcM,innerRetina,cellType)
% Initialize the spatial rf properties of a rgc mosaic given a cell type
%
% RF size scale parameters: Parasol RFs are the largest, while Midget RFs
% are about half the diameter of Parasol RFs, and SBC RFs are 1.2X the
% diameter of Parasol RFs.
%
% See Chichilnisky, E. J., and Rachel S. Kalmar. "Functional asymmetries
% in ON and OFF ganglion cells of primate retina." The Journal of
% Neuroscience 22.7 (2002).
%
% Implemented for these cell types
%
% [ON Parasol; OFF Parasol; ON Midget; OFF Midget; Small bistratified];
%
% 9/2015 JRG, ISETBIO Team
% 7/2016 JRG updated


%% Set up defaults for the sizes and weights.

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

rfSizeMult(1) = 1;      % On Parasol RF size multiplier
rfSizeMult(2) = 0.85;   % Off Parasol RF size multiplier
rfSizeMult(3) = 0.5;    % On Midget RF size multiplier
rfSizeMult(4) = 0.425;  % Off Midget RF size multiplier
rfSizeMult(5) = 1.1;    % SBC RF size multiplier

% Calculate spatial RF diameter for ON Parasol cell at a particular TEE
% See Chichilnisky, E. J., and Rachel S. Kalmar. "Functional asymmetries
% in ON and OFF ganglion cells of primate retina." The Journal of
% Neuroscience 22.7 (2002), Fig. 5, pg. 2741. 2STD fit in micrometers.
receptiveFieldDiameterParasol2STD = receptiveFieldDiameterFromTEE(innerRetina.temporalEquivEcc);

% If os is osIdentity, determine number of pixels per spatial RF diameter.
% If os is osLinear or osBioPhys, determine the number of cones per spatial
% RF diameter.
patchSizeX = innerRetina.spacing;
sensorRows = innerRetina.row;

% Needs to be a proper get from the inner retina
umPerSensorPx = patchSizeX/sensorRows;

% The spatial RF diameter in pixels (or cones) is therefore the diameter in
% microns divided by the number of microns per pixel (or cone), scaled by
% the factor determined by the type of mosaic that is being created.

% in micrometers; divide by umPerScenePx to get pixels
rgcM.rfDiameter = rfSizeMult(cellType)*(receptiveFieldDiameterParasol2STD/2); 

% Build spatial RFs of all RGCs in this mosaic
[rgcM.sRFcenter, rgcM.sRFsurround, rgcM.rfDiaMagnitude, rgcM.cellLocation, rgcM.tonicDrive] = ...
    buildSpatialRFArray(innerRetina.spacing, innerRetina.row, innerRetina.col, rgcM.rfDiameter);
% sRFcenter, sRFsurround and cellLocation are in units of the inputObj - scene pixels,
% cones or bipolar cells.
% rfDiaMagnitude is in units of micrometers.
% tonicDrive is in units of conditional intensity.


end
