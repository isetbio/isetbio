function [rfPositionsMicrons, rfPositionsDegs] = finalMRGCPositions(...
    sourceLatticeSizeDegs, eccMicrons, sizeMicrons, whichEye, ...
    MMsToDegsConversionFunction)
% Crop mRGC RF-center positions out of a retina-wide source lattice
%
% Syntax
%   [rfPositionsMicrons, rfPositionsDegs] = finalMRGCPositions(...
%        sourceLatticeSizeDegs, eccMicrons, sizeMicrons, whichEye, ...
%        MMsToDegsConversionFunction)
%
% Inputs
%   sourceLatticeSizeDegs       - Nominal generation field of view, in degrees
%   eccMicrons                  - Crop center, in retinal microns
%   sizeMicrons                 - Crop size, in retinal microns
%   whichEye                    - 'left eye' or 'right eye'
%   MMsToDegsConversionFunction - Handle converting retinal mm to degrees
%
% Returns
%   rfPositionsMicrons - Cropped positions, in retinal microns
%   rfPositionsDegs    - The same positions, in degrees
%
% Description
%   The source lattice is resolved by retinalattice.import.sourceLatticeFile,
%   which prefers a locally generated lattice and otherwise downloads the
%   published copy from the isetbio-mosaics SDR deposit into the local cache.
%
% See also
%   retinalattice, retinalattice.import.sourceLatticeFile,
%   retinalattice.import.finalConePositions

    % Load mRGC RF positions
    theMosaicFileName = retinalattice.import.sourceLatticeFile(...
        sourceLatticeSizeDegs, 'midget ganglion cells', whichEye);

    load(theMosaicFileName, 'rfPositions');

    % Reverse the polarity
    rfPositions = -rfPositions;
    
    rfPositionsMicrons = double(retinalattice.compute.croppedPositions(rfPositions, eccMicrons, sizeMicrons));
    
    % Convert positions to degs
    rfPositionsDegs = MMsToDegsConversionFunction(rfPositionsMicrons*1e-3);
end

