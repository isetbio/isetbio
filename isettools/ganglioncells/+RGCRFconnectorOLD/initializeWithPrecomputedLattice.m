function RGCRFposMicrons = initializeWithPrecomputedLattice(theInputConeMosaic)
% Initialize RGC RF positions using a precomputed lattice
%
% Syntax:
%   RGCRFposMicrons = RGCRFconnector.initializeWithPrecomputedLattice(theInputConeMosaic)
%
% Description:
%   Initialize RGC RF positions using a precomputed lattice
%
% Inputs:
%    theInputConeMosaic     - The cone mosaic that will provide input
%
% Outputs:
%    RGCRFposMicrons        - The positions (in microns) of the RGC RFs 
%
% Optional key/value pairs
%   none
%   
% History:
%   5/11/2022       NPC     Wrote it
%

    % Generate RGCRF positions that lie within the limits of the input cone mosaic
    allConePositions = theInputConeMosaic.coneRFpositionsMicrons;
    minConePosX = min(allConePositions(:,1));
    minConePosY = min(allConePositions(:,2));
    maxConePosX = max(allConePositions(:,1));
    maxConePosY = max(allConePositions(:,2));

    eccMicrons = [mean(allConePositions(:,1)) mean(allConePositions(:,2))];
    sizeMicrons = max([maxConePosX-minConePosX maxConePosY-minConePosY]);

    % Import positions from a pre-generated RGC RF lattice
    sourceLatticeSizeDegs = 58;
    customDegsToMMsConversionFunction = @(x) RGCmodels.Watson.convert.rhoDegsToMMs(x);
    
    RGCRFposMicrons = retinalattice.import.finalMRGCPositions(...
        sourceLatticeSizeDegs, ...
        eccMicrons, ...
        sizeMicrons, ...
        theInputConeMosaic.whichEye, ...
        customDegsToMMsConversionFunction); 

    idx = find(...
        (RGCRFposMicrons(:,1) >= minConePosX) & ...
        (RGCRFposMicrons(:,1) <= maxConePosX) & ...
        (RGCRFposMicrons(:,2) >= minConePosY) & ...
        (RGCRFposMicrons(:,2) <= maxConePosY));
    RGCRFposMicrons = RGCRFposMicrons(idx,:);


end