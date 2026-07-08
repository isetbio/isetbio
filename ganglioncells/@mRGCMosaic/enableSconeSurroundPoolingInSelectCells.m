function enableSconeSurroundPoolingInSelectCells(obj, theTargetRGCindices, varargin)

    % Parse input
    p = inputParser;
    p.addParameter('sConePoolingWeightsComputation', 'linearInterpolationOfLMconePoolingWeights', @ischar);
    p.addParameter('visualizeBeforeAndAfterConePoolingMaps', false, @islogical);
    p.addParameter('exportVisualizationPDFdirectory', 'mosaicComponentVisualizations', @ischar);

    % Parse input
    p.parse(varargin{:});
    sConePoolingWeightsComputation = p.Results.sConePoolingWeightsComputation;
    visualizeBeforeAndAfterConePoolingMaps = p.Results.visualizeBeforeAndAfterConePoolingMaps;
    exportVisualizationPDFdirectory = p.Results.exportVisualizationPDFdirectory;
 

    % Go through each of the RGCs for which we will enable surround S-cone pooling
    for iRGC = 1:numel(theTargetRGCindices)

        theRGCindex = theTargetRGCindices(iRGC);

        fprintf('Enabling S-cone inputs to RF surround of %d cell (%d of %d) \n', theRGCindex, iRGC ,numel(theTargetRGCindices));

        % Extract the weights and indices of cones being pooled by the RF surround
        surroundConnectivityVector = full(squeeze(obj.rgcRFsurroundConeConnectivityMatrix(:, theRGCindex)));
        surroundConeIndices = find(surroundConnectivityVector > mRGCMosaic.minSurroundWeightForInclusionInComputing);
        maxSurroundPoolingWeight = max(surroundConnectivityVector);

        % Visualize the cell's cone pooling map before enabling S-cone surround pooling
        if (visualizeBeforeAndAfterConePoolingMaps)     
            % Include surround cones whose pooling weights are >= 0.001
            minSurroundConeWeightVisualized = 0.001;
            minCenterConeWeightVisualized = mRGCMosaic.sensitivityAtPointOfOverlap;

            figPos = [1000 500];
            scaleBarDegs = 0.1;

            % Visualize cone pooling weights before enabling surround S-cone pooling
            obj.visualizeCenterSurroundConePoolingMap(theRGCindex, ...
                'minConeWeightForVisualizingRFcenterPooling', minCenterConeWeightVisualized, ...
                'minConeWeightForVisualizingRFsurroundPooling', minSurroundConeWeightVisualized, ...
                'minSurroundConeWeightRelativity', 'center', ...
                'withLineWeightingFunctions', true, ...
                'visualizedSurroundPoolingWeightMapRange', [0 maxSurroundPoolingWeight], ...
                'plotTitle', 'pre-enabling S-cone input', ...
                'scaleBarDegs', scaleBarDegs, ...
                'doNotLabelScaleBar', true, ...
                'figNo', 100, ...
                'figPos', figPos, ...
                'exportToFigurePDFsDirWithPDFFileName', sprintf('RFmap%d_PreSconeInput.pdf', theRGCindex), ...
                'pdfExportSubDir', exportVisualizationPDFdirectory);
        end

        switch (sConePoolingWeightsComputation)

            case 'linearInterpolationOfLMconePoolingWeights'
                surroundIntegratedSensitivityBeforeEnablingSconeInputs = sum(surroundConnectivityVector(surroundConeIndices));
                
        
                surroundPooledLMconesXcoords = obj.inputConeMosaic.coneRFpositionsDegs(surroundConeIndices,1);
                surroundPooledLMconesYcoords = obj.inputConeMosaic.coneRFpositionsDegs(surroundConeIndices,2);
                surroundPooledLMconesWeights = surroundConnectivityVector(surroundConeIndices);
        
                F = scatteredInterpolant(surroundPooledLMconesXcoords(:),surroundPooledLMconesYcoords(:),surroundPooledLMconesWeights(:));
        
                % Compute the radius of surround pooling for this cell
                surroundRadiusDegs = max(sqrt(sum(bsxfun(@minus, obj.inputConeMosaic.coneRFpositionsDegs(surroundConeIndices,:), obj.rgcRFpositionsDegs(theRGCindex,:)).^2,2)));
                
                % Compute distances of all S-cones in the input cone mosaic to the RF center
                sConeDistancesToRGCRFcenter = sqrt(sum(bsxfun(@minus, obj.inputConeMosaic.coneRFpositionsDegs(obj.inputConeMosaic.sConeIndices,:), obj.rgcRFpositionsDegs(theRGCindex,:)).^2,2));
        
                % Determine the indices of S-cones whose distance to the RF center is
                % within the radius of surround pooling for this RGC
                sConeIndicesConsidered = obj.inputConeMosaic.sConeIndices(find(sConeDistancesToRGCRFcenter < surroundRadiusDegs));
        
 
                if (numel(sConeIndicesConsidered) == 0)
                    fprintf('No S-cones in the input cone mosaic surrounding mRGC %d\n', theRGCindex);
                else
                
                    % Interpolate S-cone weights for all sConeIndicesConsidered
                    theSconeWeights = F(obj.inputConeMosaic.coneRFpositionsDegs(sConeIndicesConsidered,1), obj.inputConeMosaic.coneRFpositionsDegs(sConeIndicesConsidered,2));
            
                    % Add S-cone weights
                    obj.rgcRFsurroundConeConnectivityMatrix(sConeIndicesConsidered, theRGCindex) = theSconeWeights;
            
                    surroundConnectivityVector = full(squeeze(obj.rgcRFsurroundConeConnectivityMatrix(:, theRGCindex)));
                    surroundConeIndices = find(surroundConnectivityVector > mRGCMosaic.minSurroundWeightForInclusionInComputing);
                    surroundIntegratedSensitivityAfterEnablingSconeInputs = sum(surroundConnectivityVector(surroundConeIndices));
            
                    % Factor to match total surround weight
                    surroundSensitivityMatchingFactor = surroundIntegratedSensitivityBeforeEnablingSconeInputs/surroundIntegratedSensitivityAfterEnablingSconeInputs;
            
                    % Scale surround weights according to surroundSensitivityMatchingFactor
                    obj.rgcRFsurroundConeConnectivityMatrix(:, theRGCindex) = ...
                        obj.rgcRFsurroundConeConnectivityMatrix(:, theRGCindex) * surroundSensitivityMatchingFactor;
               end

            otherwise
                error('Unkown sConePoolingWeightsComputation method: ''%s''.', sConePoolingWeightsComputation);

        end  %  switch (sConePoolingWeightsComputation)

        % Visualize the cell's cone pooling map after enabling S-cone surround pooling
        if (visualizeBeforeAndAfterConePoolingMaps)
             % Visualize cone pooling weights after enabling surround S-cone pooling
             obj.visualizeCenterSurroundConePoolingMap(theRGCindex, ...
                'minConeWeightForVisualizingRFcenterPooling', minCenterConeWeightVisualized, ...
                'minConeWeightForVisualizingRFsurroundPooling', minSurroundConeWeightVisualized, ...
                'minSurroundConeWeightRelativity', 'center', ...
                'withLineWeightingFunctions', true, ...
                'visualizedSurroundPoolingWeightMapRange', [0 maxSurroundPoolingWeight], ...
                'plotTitle', 'post-enabling S-cone input', ...
                'scaleBarDegs', scaleBarDegs, ...
                'doNotLabelScaleBar', true, ...
                'figNo', 101, ...
                'figPos', figPos + [200 0], ...
                'exportToFigurePDFsDirWithPDFFileName', sprintf('RFmap%d_PostSconeInput.pdf', theRGCindex), ...
                'pdfExportSubDir', exportVisualizationPDFdirectory);
        end

    end % iRGC
end
