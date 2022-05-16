function [RGCRFinputs, RGCRFweights] = wireInputConeMosaicToRGCRFs(...
        theInputConeMosaic, RGCRFposMicrons, ...
        wiringParams,visualizationParams)
% Wire the input cone mosaic to the RF centers of a set of RGCs
%
% Syntax:
%   [RGCRFinputs, RGCRFweights] = RGCRFconnector.wireInputConeMosaicToRGCRFs(...
%        theInputConeMosaic, RGCRFposMicrons, ...
%        chromaticSpatialVarianceTradeoff, ...
%        sequentialConeTypeWiring, ...
%        maxNearbyRGCsNum, ...
%        generateProgressionVideo)
%
% Description:
%   Wire the input cone mosaic to the RF centers of a set of RGCs
%
% Inputs:
%    theInputConeMosaic                 - the input cone mosaic
%    RGCRFposMicrons                    - [M x 2] matrix of (x,y) positions (in microns) of M target RGC RF centers
%    wiringParams                       - Struct with the following wiring params
%       chromaticSpatialVarianceTradeoff   - Chromatic-SpatialVariance tradefoff
%       sequentialConeTypeWiring           - Logical. Whether to wire cone types sequentially or not
%       maxNearbyRGCsNum                   - Scalar. Max number of nearby RGCs to look for assignment, if a cone cannot
%                                         be assigned to its closest RGC because that RGC already has one input cone
%    visualizationParams                 - Struct with the following visualization params
%       generateProgressionVideo           - Boolean. Whether to generate a video showning the different stages of the wiring
%
% Outputs:
%    RGCRFinputs                - Cell array with indices of the input cones, one cell per each target RGC RF center
%    RGCRFweights               - Cell array with weights of the input cones, one cell per each target RGC RF center   
%
% Optional key/value pairs
%   none
%   
% History:
%   5/11/2022       NPC     Wrote it
%

    

    % Step1. Connect cones to the nearestRGC that minimizes the cost
    % according to the passed chromaticSpatialVarianceTradeoff
    [RGCRFinputs, RGCRFweights, availableZeroInputRGCsNum] = ...
        RGCRFconnector.connectEachConeToItsNearestRGC(...
            theInputConeMosaic, RGCRFposMicrons, ...
            wiringParams, visualizationParams);

    % Step2. 

end