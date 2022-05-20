function ratios = localConeToRGCDensityRatiosAtScatteredPositions(...
        localConeToRGCDensityRatioStruct, RGCRFposMicrons)
% Compute the local cone/RGC density ratios at some RGC position(s)
%
% Syntax:
%   ratios = RGCRFconnector.localConeToRGCDensityRatiosAtScatteredPositions(...
%        localConeToRGCDensityRatioStruct, RGCRFpos);
%
% Description:
%   Compute the local cone/RGC density ratios at some RGC position(s),
%   specified in microns
%
% Inputs:
%    localConeToRGCDensityRatioStruct   - Struct with info to compute the cone/RGC density ratio at any (x,y) position
%    RGCRFposMicrons                    - [M x 2] matrix of (x,y) positions (in microns) of M RGC RF centers
%
% Outputs:
%    ratios                             - Vector with computed cone/RGC density ratios
% 
% Optional key/value pairs
%   none
%   
% History:
%   5/11/2022       NPC     Wrote it
%

    X = RGCRFposMicrons(:,1);
    Y = RGCRFposMicrons(:,2);
    coneDensities = localConeToRGCDensityRatioStruct.coneDensityFunctionHandle(X,Y);
    rgcDensities = localConeToRGCDensityRatioStruct.RGCDensityFunctionHandle(X,Y);
    ratios = coneDensities ./ rgcDensities;

end