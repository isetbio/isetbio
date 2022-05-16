function [cost, rfCentroid, spatialVarianceCost, chromaticVarianceCost] = costToMaintainInputCones(w, ...
    inputConePositions, inputConeSpacings, ...
    inputConeTypes, inputConeWeights)
% Compute the cost of maintaining a set of cone inputs 
%
% Syntax:
%   [cost, rfCentroid] = RGCRFconnector.costToMaintainInputCones(w, ...
%    inputConePositions, inputConeSpacings, ...
%    inputConeTypes, inputConeWeights)
%
% Description:
%     Compute the cost of maintaining a set of cone inputs. The cost is a function
%     of spatial and/or chromaticthere. The balance between the two is 
%     controlled by the value of w (chromaticSpatialVarianceTradeoff), as follows:
%
% w      COST is ...                                                    Operation 
% 1      proportional to spatial variance only                          Minimize spatial variance, ignoring homogeneity of input cone types
% 0      proportional to chromatic variance only, e.g.:
%        lConeSignal/mConeSignal, if lConeSignal < mConeSignal, or      Minimize chromatic variance, ignoring spatial position of input cones
%        mConeSignal/lConeSignal, if mConeSignal < lConeSignal
% 0<w<1  w * spatialVariance + (1-w)*chromatic variance                 Minimize the combined chromatic+spatial variance
%
%
% Inputs:
%    w                      - Chromatic-SpatialVariance tradefoff
%    inputConePositions     - Positions of the input cones
%    inputConeSpacings      - Spacings of the input cones
%    inputConeTypes         - Types of the input cones
%    inputConeWeights       - Weights of the input cones
%
% Outputs:
%    cost                   - The total cost of maintaining the cone inputs
%    rfCentroid             - The geometric centroid of the RGCRF based on
%                             the weighted spatial positions of the input cones
%    spatialVarianceCost    - The spatial variance cost component
%    chromaticVarianceCost  - The chromatic variance cost component
%
% History:
%   5/11/2022       NPC     Wrote it
%

    assert(size(inputConePositions,1)>=1, 'RGC must have at least one input cone.');
    assert(size(inputConePositions,2)==2, 'input cone positions must be [nx2]');
    assert(((w>=0)&&(w<=1)), 'w must be between 0 and 1')

    if (size(inputConePositions,1) == 1)
        cost = 0;
        rfCentroid = inputConePositions;
        spatialVarianceCost = 0;
        chromaticVarianceCost = 0;
        return;
    end

    if (isempty(inputConeWeights))
        inputConeWeights = ones(size(inputConePositions,1), 1);
    else
        inputConeWeights = reshape(inputConeWeights, [numel(inputConeWeights) 1]);
    end

    % SPATIAL VARIANCE COST
    measure = 'max distance';
    measure = 'weighted variance';
    switch (measure)
        case 'max distance'
            allDistances = pdist2(inputConePositions, inputConePositions);
            spatialVarianceCost = max(allDistances(:))/mean(inputConeSpacings);
            [~, rfCentroid] = var(inputConePositions,inputConeWeights,1);

        case 'weighted variance'
            % Compute weighted variance
            [spatialVarianceCost, rfCentroid] = var(inputConePositions,inputConeWeights,1);
            spatialVarianceCost = sqrt(mean(spatialVarianceCost))/mean(inputConeSpacings);

        otherwise
            error('unknown measure')
    end


    % CHROMATIC VARIANCE COST
    % non-dominant cone-type / dominant cone-type
    lConeIndices = find(inputConeTypes == cMosaic.LCONE_ID);
    mConeIndices = find(inputConeTypes == cMosaic.MCONE_ID);

    lConeSignal = sum(inputConeWeights(lConeIndices));
    mConeSignal = sum(inputConeWeights(mConeIndices));
    totalLMConeSignal = lConeSignal + mConeSignal;
    if (lConeSignal <= mConeSignal)
        chromaticVarianceCost = lConeSignal/totalLMConeSignal;
    else
        chromaticVarianceCost = mConeSignal/totalLMConeSignal;
    end

    % Total cost
    cost = w * spatialVarianceCost + (1-w)*chromaticVarianceCost;
end

