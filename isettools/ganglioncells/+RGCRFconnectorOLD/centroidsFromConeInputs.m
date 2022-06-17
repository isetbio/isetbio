function [centroids, coneInputsNum] = centroidsFromConeInputs(RGCRFinputs, RGCRFweights, coneRFpos)
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
%    coneInputsNum      - [N x1] vector of # of cone inputs for each RGC
%
% Optional key/value pairs
%   none
%   
% History:
%   5/11/2022       NPC     Wrote it
%

    rgcsNum = numel(RGCRFinputs);
    centroids = zeros(rgcsNum,2);
    coneInputsNum = zeros(numel(RGCRFinputs),1);
    

    parfor iRGC = 1:rgcsNum
        coneInputsNum(iRGC) = numel(RGCRFinputs{iRGC});
        inputConePositions = coneRFpos(RGCRFinputs{iRGC},:);
        if (isempty(RGCRFweights))
            inputConeWeights = ones(1,numel(RGCRFinputs{iRGC}));
        else
            inputConeWeights = RGCRFweights{iRGC};
            if (isempty(inputConeWeights))
                inputConeWeights = ones(1,numel(RGCRFinputs{iRGC}));
            end
        end

        [~, centroids(iRGC,:)] = var(inputConePositions,inputConeWeights,1);
    end    

end