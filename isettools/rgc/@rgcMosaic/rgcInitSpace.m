function rgcM = rgcInitSpace(rgcM,innerRetina,cellType)
% Initialize the spatial rf properties of a rgc mosaic for a cell type
%
%    rgcM = rgcInitSpace(rgcM,innerRetina,cellType)
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
%   ir = irCreate(bipolar(coneMosaic));
%   ir.mos
% See also:
%
% JRG/BW ISETBIO Team 2015

%% All three arguments are required

p = inputParser;
p.addRequired('rgcM');
p.addRequired('innerRetina');

vFunc = @(x)(ismember(ieParamFormat(x),{'onparasol', 'offparasol', 'onmidget', 'offmidget', 'smallbistratified'}));
p.addRequired('cellType',vFunc);

p.parse(rgcM,innerRetina,cellType);

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

% Calculate spatial RF diameter for ON Parasol cell at a particular TEE
% See Chichilnisky, E. J., and Rachel S. Kalmar. "Functional asymmetries
% in ON and OFF ganglion cells of primate retina." The Journal of
% Neuroscience 22.7 (2002), Fig. 5, pg. 2741. 2STD fit in micrometers.
receptiveFieldDiameterParasol2STD = ...
    receptiveFieldDiameterFromTEE(innerRetina.temporalEquivEcc);

% The spatial RF diameter in pixels (or cones) is therefore the diameter in
% microns divided by the number of microns per pixel (or cone), scaled by
% the factor determined by the type of mosaic that is being created.

% in micrometers; divide by umPerScenePx to get pixels
rgcM.rfDiameter = rfSizeMult*(receptiveFieldDiameterParasol2STD/2); 

% Build spatial RFs of all RGCs in this mosaic
[rgcM.sRFcenter, rgcM.sRFsurround, rgcM.rfDiaMagnitude, rgcM.cellLocation, rgcM.tonicDrive, rgcM.ellipseMatrix] = ...
    buildSpatialRFArray(innerRetina.size, innerRetina.row, innerRetina.col, rgcM.rfDiameter);

% sRFcenter, sRFsurround and cellLocation are in units of the inputObj - scene pixels,
% cones or bipolar cells.
% rfDiaMagnitude is in units of micrometers.
% tonicDrive is in units of conditional intensity.

end
