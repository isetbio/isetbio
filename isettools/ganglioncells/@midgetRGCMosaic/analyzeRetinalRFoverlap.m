function [NNNDs, NNNDtuplets, RGCdistances, distancesFromMosaicCenterDegs, targetRGCindices] = analyzeRetinalRFoverlap(obj, varargin)
    p = inputParser;
    p.addParameter('spatialSupportSamplesNum', 256, @isnumeric);
    p.addParameter('neighborsNumForAnalysis', 3, @isnumeric);
    p.parse(varargin{:});
    spatialSupportSamplesNum = p.Results.spatialSupportSamplesNum;
    neighborsNumForAnalysis = p.Results.neighborsNumForAnalysis;

    % Compute additional margin for the spatial support
    mRGCmosaicCenterDegs = mean(obj.rgcRFpositionsDegs,1);
    marginDegs = min([0.5 0.2*min(obj.sizeDegs)]);

     % Compute the retinal RFcenter maps
    theRetinalRFcenterMaps = obj.computeRetinalRFcenterMaps(marginDegs, spatialSupportSamplesNum);

    rgcsNum = size(obj.rgcRFpositionsDegs,1);

    % Fitting the discrete RF center cone map
    fitTheDiscreteRFcenterMap = true;

    % Fit an ellipsoidal Gaussian and extract the sigma to all the RGCs
    sigmaDegs = zeros(1, rgcsNum);
    parfor iRGC = 1:rgcsNum
        % Retrieve the computed retinal center RF map
        s = theRetinalRFcenterMaps{iRGC};

        if (fitTheDiscreteRFcenterMap)
            % Fit the discrete center RF map with an ellipsoidal Gaussian
            theFittedGaussian = RetinaToVisualFieldTransformer.fitScatterGaussianEllipsoid(...
                s.spatialSupportDegsX, s.spatialSupportDegsY, s.centerRF,...
                s.inputConeWeights, obj.inputConeMosaic.coneRFpositionsDegs(s.inputConeIndices,:), ...
                'flatTopGaussian', ~true, ...
                'forcedOrientationDegs', [], ...
                'rangeForEllipseRcYRcXratio', [1/4 4], ...
                'forcedCentroidXYpos', obj.rgcRFpositionsDegs(iRGC,:), ...
                'globalSearch', true, ...
                'multiStartsNum', 8);
        else
            % Fit the continuous center RF map with an ellipsoidal Gaussian
            theFittedGaussian = RetinaToVisualFieldTransformer.fitGaussianEllipsoid(...
                s.spatialSupportDegsX, s.spatialSupportDegsY, s.centerRF, ...
                'flatTopGaussian', ~true, ...
                'forcedOrientationDegs', [], ...
                'rangeForEllipseRcYRcXratio', [1/1.4 1.4], ...
                'forcedCentroidXYpos', obj.rgcRFpositionsDegs(iRGC,:), ...
                'globalSearch', true, ...
                'multiStartsNum', 8);
        end

        % Rc = sqrt(RcX * RcY)
        characteristicRadiusDegs = sqrt(prod(theFittedGaussian.characteristicRadii));

        % sigma = Rc / sqrt(2.0)
        sigmaDegs(iRGC) = characteristicRadiusDegs/sqrt(2.0);
    end

    % Compute the NormalizedNearestNeighborDistance as per Gauthier, Chichilinsky et al (2009)
    NNNDtuplets = nan(rgcsNum,neighborsNumForAnalysis);
    RGCdistances = nan(rgcsNum,neighborsNumForAnalysis);

    for targetRGCindex = 1:rgcsNum
         targetPosDegs = obj.rgcRFpositionsDegs(targetRGCindex,:);
         targetSigmaDegs = sigmaDegs(targetRGCindex);

         % Find the neighboring RGCs
         [~, idx] = MosaicConnector.pdist2(obj.rgcRFpositionsDegs, targetPosDegs, ...
                    'smallest', 1+neighborsNumForAnalysis);
         theNeighboringRGCindices = idx(2:end);

         for iNeighbor = 1:numel(theNeighboringRGCindices)
             theNeighborRGCindex = theNeighboringRGCindices(iNeighbor);
             neighborPosDegs = obj.rgcRFpositionsDegs(theNeighborRGCindex ,:);
             neighborSigmaDegs = sigmaDegs(theNeighborRGCindex);
             centerSeparation = sqrt(sum((targetPosDegs-neighborPosDegs).^2,2));
             NNNDtuplets(targetRGCindex, iNeighbor) = 2*centerSeparation/(targetSigmaDegs+neighborSigmaDegs);
             RGCdistances(targetRGCindex, iNeighbor) = centerSeparation;
         end

    end

    % Sort RGCs according to their eccentricity
    ecc = sqrt(sum((bsxfun(@minus, obj.rgcRFpositionsDegs, mRGCmosaicCenterDegs)).^2,2));
    targetRGCindices = find(ecc < 0.8*max(obj.sizeDegs)*0.5);
    distancesFromMosaicCenterDegs = ecc(targetRGCindices);
    NNNDtuplets = NNNDtuplets(targetRGCindices, :);
    RGCdistances = RGCdistances(targetRGCindices, :);
    fprintf('Analyzing the NNNDs of the central %d RGCs out of a total of %d RGCs\n', numel(targetRGCindices), rgcsNum)


    NNNDs = [];
    for iRGC = 1:size(NNNDtuplets,1)
        % Find the valid neigbor indices
        iNeighborIndices = find(~isnan(squeeze(NNNDtuplets(iRGC,:))));
        % Accumulate the NNNDs
        NNNDs(numel(NNNDs) + (1:numel(iNeighborIndices))) = NNNDtuplets(iRGC,iNeighborIndices);
    end
    
end