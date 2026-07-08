function cost = spatialVarianceCost(spatialVarianceMetric, inputWeights, inputPositions,  destinationRFspacing)

    switch (spatialVarianceMetric)
            case 'maximal interinput distance'
                spatialVarianceCost = MosaicConnector.maximalInterInputDistance(inputPositions);
            case 'spatial variance'
                varianceXY = var(inputPositions,inputWeights,1);
                spatialVarianceXY = varianceXY(:);
                spatialVarianceCost = sqrt(spatialVarianceXY(1)+spatialVarianceXY(2));
            otherwise
                error('Unknown spatialVarianceMetric: ''%s''.', obj.wiringParams.spatialVarianceMetric);
    end

    cost = spatialVarianceCost / destinationRFspacing;
end
