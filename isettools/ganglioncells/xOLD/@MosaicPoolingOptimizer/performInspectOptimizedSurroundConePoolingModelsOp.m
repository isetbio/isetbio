function performInspectOptimizedSurroundConePoolingModelsOp(mosaicParams, varargin)

   % Parse optional input
    p = inputParser;
    p.addParameter('tickSeparationArcMin', [], @(x)(isempty(x)||isscalar(x)));
    p.addParameter('normalizedPeakSurroundSensitivity', 0.4, @isscalar);
    p.addParameter('visualizedSpatialFrequencyRange', [], @(x)(isempty(x)||(numel(x)==2)));
    p.addParameter('gridlessLineWeightingFuncions', false, @islogical);
    p.parse(varargin{:});

    tickSeparationArcMin = p.Results.tickSeparationArcMin;
    normalizedPeakSurroundSensitivity = p.Results.normalizedPeakSurroundSensitivity;
    visualizedSpatialFrequencyRange = p.Results.visualizedSpatialFrequencyRange;
    gridlessLineWeightingFuncions = p.Results.gridlessLineWeightingFuncions;

    % Generate the mosaic filename
    [mosaicFileName, resourcesDirectory] = ...
        MosaicPoolingOptimizer.resourceFileNameAndPath('mosaic', ...
            'mosaicParams', mosaicParams);

    % Load the generated center-only connected mRGCmosaic
    load(fullfile(resourcesDirectory,mosaicFileName), 'theMidgetRGCMosaic');


    % Ask the user which optics were used for computing the input cone
    % mosaic STF responses
    opticsParams = MosaicPoolingOptimizer.chooseOpticsForInputConeMosaicSTFresponses(mosaicParams);

    % Ask the user which H1 cell index was used to optimize the RF
    % surround pooling model models
    [~, gridSamplingScheme,  optimizedRGCpoolingObjectsFileName] = ...
            MosaicPoolingOptimizer.chooseRFmodelForSurroundConePoolingOptimization(mosaicParams, opticsParams);


    % Instantiate a MosaicPoolingOptimizer object with the center-connected
    % mRGC mosaic
    theMosaicPoolingOptimizerOBJ = MosaicPoolingOptimizer(...
         theMidgetRGCMosaic, ...
         'samplingScheme', gridSamplingScheme, ...
         'generateSamplingGrids', true, ...
         'visualizeSamplingGrids', false);

    % Visualize the grid nodes
    visualizeNodes(theMosaicPoolingOptimizerOBJ, optimizedRGCpoolingObjectsFileName);


    queryString = sprintf('\nEnter grid node to inspect [%d-%d]. Alternatively, hit enter to inspect multiple nodes: ', ...
        1, theMosaicPoolingOptimizerOBJ.gridNodesNum);
    gridNodesToInspect = input(queryString, 's');
    if (isempty(gridNodesToInspect))
        gridNodesToInspect = MosaicPoolingOptimizer.gridNodesToOptimize();
        if (isempty(gridNodesToInspect))
             gridNodesToInspect = 1:theMosaicPoolingOptimizerOBJ.gridNodesNum;
        end
    else
       gridNodesToInspect = str2num(gridNodesToInspect);
    end


    for iNode = 1:numel(gridNodesToInspect)
       if (iscell(gridNodesToInspect))
            gridNodeIndex = gridNodesToInspect{iNode}.number;
       else
            gridNodeIndex = gridNodesToInspect(iNode);
       end

       theMosaicPoolingOptimizerOBJ.inspect(...
           gridNodeIndex, ...
           opticsParams, ...
           optimizedRGCpoolingObjectsFileName, ...
           'tickSeparationArcMin', tickSeparationArcMin, ...
           'normalizedPeakSurroundSensitivity', normalizedPeakSurroundSensitivity, ...
           'visualizedSpatialFrequencyRange', visualizedSpatialFrequencyRange, ...
           'gridlessLineWeightingFuncions', gridlessLineWeightingFuncions ...
           );

       fprintf('Exported files for node %d (%d of %d nodes)\n', gridNodeIndex, iNode , numel(gridNodesToInspect));
       fprintf('Hit enter to continue ...');
       pause
    end

end


function  visualizeNodes(theMosaicPoolingOptimizerOBJ, optimizedRGCpoolingObjectsFileName)

    fprintf('Loading coneRFcomputeStructs for all mosaic nodes to visualize their positions\n');

    allGridNodes = 1:theMosaicPoolingOptimizerOBJ.gridNodesNum;
    gridPositions = nan(2,max(theMosaicPoolingOptimizerOBJ.conesNumPooledByTheRFcenterGrid), numel(allGridNodes),2);

    for gridNodeIndex = 1:numel(allGridNodes)
        optimizedRGCpoolingObjectsFileNameForThisNode = ...
               strrep(optimizedRGCpoolingObjectsFileName, '.mat', sprintf('_ForGridNode_%d.mat', gridNodeIndex));

        load(optimizedRGCpoolingObjectsFileNameForThisNode, ...
            'theLconeRFcomputeStruct', ...
            'theMconeRFcomputeStruct');

        % Check that the number of center cones in the grid location match
        % those in the compute structs
        conesNum = theMosaicPoolingOptimizerOBJ.conesNumPooledByTheRFcenterGrid(gridNodeIndex);
        assert(conesNum == numel(theLconeRFcomputeStruct.modelConstants.indicesOfCenterCones), 'Inconsistency of cones num for LconeRFcompute struct');
        assert(conesNum == numel(theMconeRFcomputeStruct.modelConstants.indicesOfCenterCones), 'Inconsistency of cones num for MconeRFcompute struct');

        % Check that the majority type of center cones in the grid location match
        % that in the compute structs
        coneTypes = theLconeRFcomputeStruct.modelConstants.inputConeMosaicConeTypes(theLconeRFcomputeStruct.modelConstants.indicesOfCenterCones);
        LconesNum = numel(find((coneTypes == cMosaic.LCONE_ID) == 1));
        MconesNum = numel(find((coneTypes == cMosaic.MCONE_ID) == 1));
        assert(LconesNum >= MconesNum, 'Inconsistency of majority cones for LconeRFcompute struct');
        coneTypes = theMconeRFcomputeStruct.modelConstants.inputConeMosaicConeTypes(theMconeRFcomputeStruct.modelConstants.indicesOfCenterCones);
        LconesNum = numel(find((coneTypes == cMosaic.LCONE_ID) == 1));
        MconesNum = numel(find((coneTypes == cMosaic.MCONE_ID) == 1));
        assert(MconesNum >= LconesNum, 'Inconsistency of majority cones for MconeRFcompute struct');

        % Get grid location 
        LconeRGCindex = theMosaicPoolingOptimizerOBJ.targetRGCindicesWithLconeMajorityCenter(gridNodeIndex);
        MconeRGCindex = theMosaicPoolingOptimizerOBJ.targetRGCindicesWithMconeMajorityCenter(gridNodeIndex);
        LcenterRFposition = theMosaicPoolingOptimizerOBJ.theRGCMosaic.rgcRFpositionsDegs(LconeRGCindex,:);
        McenterRFposition = theMosaicPoolingOptimizerOBJ.theRGCMosaic.rgcRFpositionsDegs(MconeRGCindex,:);

        gridPositions(1, conesNum, gridNodeIndex,:) = LcenterRFposition;
        gridPositions(2, conesNum, gridNodeIndex,:) = McenterRFposition;
    end

    theConesNumCases = unique(theMosaicPoolingOptimizerOBJ.conesNumPooledByTheRFcenterGrid);
    hFig = figure(100); clf;
    set(hFig, 'Color', [0 0 0]);
    xPos = gridPositions(:, :, :,1);
    yPos = gridPositions(:, :, :,2);
    xLims = [min(xPos(:))-0.5 max(xPos(:))+0.5];
    yLims = [min(yPos(:))-0.5 max(yPos(:))+0.5];

    for iConesNum = 1:numel(theConesNumCases)
        centerConesNum = theMosaicPoolingOptimizerOBJ.conesNumPooledByTheRFcenterGrid(iConesNum);
        ax = subplot(3,2, iConesNum);
        hold on;
        for gridNodeIndex = 1:size(gridPositions,3)
            xyPosL = squeeze(gridPositions(1, centerConesNum, gridNodeIndex,:));
            if (~isnan(xyPosL))
                plot(ax, xyPosL(1), xyPosL(2), 'ro');
            else
                continue;
            end

            xyPosM = squeeze(gridPositions(2, centerConesNum, gridNodeIndex,:));
            if (~isnan(xyPosM))
                plot(ax, xyPosM(1), xyPosM(2), 'go');
            else
                continue;
            end
            
            xyPos = mean([xyPosL xyPosM],2);
            plot(ax, [xyPosL(1) xyPos(1)], [xyPosL(2) xyPos(2)], 'r-');
            plot(ax, [xyPosM(1) xyPos(1)], [xyPosM(2) xyPos(2)], 'g-');
            text(ax,xyPos(1), xyPos(2), sprintf(' %d', gridNodeIndex), 'Color', 'w', 'FontSize', 16);
        end
        set(ax, 'XLim', xLims, 'YLim', yLims, 'Color', [0 0 0], 'XColor', [0.5 0.5 0.5], 'YColor', [0.5 0.5 0.5]);
        box(ax, 'on');
        title(sprintf('%d center cones', centerConesNum), 'Color', [1 1 1], 'FontSize', 20);
        drawnow;
    end
end
