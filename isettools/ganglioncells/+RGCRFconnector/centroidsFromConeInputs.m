function centroids = centroidsFromConeInputs(RGCRFinputs, RGCRFweights, coneRFpos)
% Compute centroids of RGC RFs fron their cone inputs
%
% Syntax:
%   centroids = RGCRFconnector.centroidsFromConeInputs(RGCRFinputs, RGCRFweights, coneRFpos)
%
% Description:
%   Compute centroids of RGC RFs fron their cone inputs
%
% Inputs:
%    RGCRFinputs        - Cell array with indices of the input cones, one cell per each target RGC RF center
%    RGCRFweights       - Cell array with weights of the input cones, one cell per each target RGC RF center 
%
% Outputs:
%    centroids          - [N x 2] matrix of (x,y) centroids of the RGC RF centers 
%
% Optional key/value pairs
%   none
%   
% History:
%   5/11/2022       NPC     Wrote it
%

    centroids = zeros(numel(RGCRFweights),2);
    parfor iRGC = 1:numel(RGCRFweights)
        inputConePositions = coneRFpos(RGCRFinputs{iRGC},:);
        inputConeWeights = RGCRFweights{iRGC};
        [~, centroids(iRGC,:)] = var(inputConePositions,inputConeWeights,1);
    end    

end