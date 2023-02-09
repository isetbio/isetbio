function RFcomputeStruct = retinalConePoolingParamsForTargetVisualRF(obj, ...
                centerConeType, initialRetinalConePoolingParamsStruct)

  
    RFcomputeStruct = [];

    indicesOfConesPooledByTheRFcenter = obj.targetVisualRFDoGparams.indicesOfConesPooledByTheRFcenter;
    weightsOfConesPooledByTheRFcenter = obj.targetVisualRFDoGparams.weightsOfConesPooledByTheRFcenter;

    switch (centerConeType)
        case cMosaic.LCONE_ID
            theRFCenterConeMajorityPSF  = obj.spectrallyWeightedPSFData.LconeWeighted;

            spatialPositionDegs = mean(obj.coneMosaic.coneRFpositionsDegs(indicesOfConesPooledByTheRFcenter,:),1);
            figureName = sprintf('RF center with %d L-cone(s) at (%2.2f,%2.2f) degs', ...
                numel(indicesOfConesPooledByTheRFcenter), spatialPositionDegs(1), spatialPositionDegs(2));
            summaryFigNo = 2000 + numel(weightsOfConesPooledByTheRFcenter)*10+1;
            
        case cMosaic.MCONE_ID
            theRFCenterConeMajorityPSF  = obj.spectrallyWeightedPSFData.MconeWeighted;

            spatialPositionDegs = mean(obj.coneMosaic.coneRFpositionsDegs(indicesOfConesPooledByTheRFcenter,:),1);
            figureName = sprintf('RF center with %d M-cone(s) at (%2.2f,%2.2f) degs', ...
                numel(indicesOfConesPooledByTheRFcenter), spatialPositionDegs(1), spatialPositionDegs(2));
            summaryFigNo = 2000 + numel(weightsOfConesPooledByTheRFcenter)*10+2;

        otherwise
            error('Not L or M cone: %d', centerConeType);
    end

end