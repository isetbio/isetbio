classdef RTVFmultifocal < handle

    % Public properties (read-only)
    properties (SetAccess = private)
        % The input midget RGC mosaic
        theRGCMosaic;

        % The input params
        mosaicCenterParams;
        opticsParams;
        rfModelParams; 

        % The compute method and its resources
        stfComputeMethod;
        stfComputeMethodResources;

        % The nominal spatial sampling grid
        nominalSpatialSamplingGrid;

        % The full sampling grids
        opticalPositionGrid;
        conesNumPooledByTheRFcenterGrid;
        visualSTFSurroundToCenterRcRatioGrid;
        visualSTFSurroundToCenterIntegratedSensitivityRatioGrid;

        % The target RGC index for each RTVFobj
        targetRGCindex;

        % The list with all computed RTVF objects
        RTVFTobjList;

        % Multistarts num
        multiStartsNumRetinalPooling;
    end % Read-only



    % Public methods (class interface)
    methods
        % Constructor
        function obj = RTVFmultifocal(theRGCmosaic, ...
            mosaicCenterParams, opticsParams, rfModelParams, ...
            samplingScheme, stfComputeMethod, varargin)

            obj.theRGCMosaic = theRGCmosaic;
            obj.mosaicCenterParams = mosaicCenterParams;
            obj.opticsParams = opticsParams;
            obj.rfModelParams = rfModelParams;
        
            assert(ismember(stfComputeMethod, RTVF.validSTFcomputeMethods), ...
                sprintf('STF computeMethod ''%s'' is not valid', stfComputeMethod));
            obj.stfComputeMethod = stfComputeMethod;

            if (strcmp(obj.stfComputeMethod, RTVF.directSTFcomputeMethod))
                % Ensure that the file with the mosaic responses exists
                midgetRGCMosaicInspector.say('Select file with cone mosaic STF responses');

                dropboxDir = midgetRGCMosaicInspector.localDropboxPath();
                [file,path] = uigetfile(fullfile(dropboxDir, '*.mat'), ...
                            'Select a file');
    
                if (file ~= 0)
                    coneMosaicResponsesFileName = fullfile(path,file);
                    obj.loadComputeMethodResources(coneMosaicResponsesFileName);
                end
            end

            % Parse optional input
            p = inputParser;
            p.addParameter('multiStartsNumRetinalPooling', [], @(x)(isempty(x)||isscalar(x)));
            p.addParameter('minSpatialSamplingDegs', 0.25, @isnumeric);
            p.parse(varargin{:});
            obj.multiStartsNumRetinalPooling = p.Results.multiStartsNumRetinalPooling;
            
            % Generate the nominal spatial sampling grid
            obj.generateNominalSpatialSamplingGrid(...
                samplingScheme, true);

            % Generate the full sampling grids
            obj.generateSamplingGrids(p.Results.minSpatialSamplingDegs, true);

        end % Constructor

        % Compute method which fits an RTVF to each sampling position for
        % each each case of # of input cones in the RF center
        compute(obj, initialGridRetinalConePoolingParamsStruct, ...
            RTVobjIndicesToBeComputed, computeLconeCenterComputeStruct, ...
            computeMconeCenterComputeStruct, exportsDirectory);

    end % Public methods

    % Private methods
    methods (Access=private)
        % Method to generate the nominal spatial sampling grid
        generateNominalSpatialSamplingGrid(obj, ...
            samplingScheme,visualizeSamplingGrid);

        % Method to generate the full sampling grids
        generateSamplingGrids(obj, minSpatialSamplingDegs, visualizeSpatialSamplingGrids);

        % Method to plot various spatial sampling grids
        plotSpatialSamplingGrid(obj, ax, spatialSamplingGrid, plotTitle);
    end % Private methods


     % Class methods
    methods (Static)
        % Method to generate the spatial sampling grid of the multi-focal RTVF
        gridCoords  = spatialSamplingGridCoords(eccentricityDegs, sizeDegs, ...
            rgcRFpositionsDegs, gridHalfSamplesNum, samplingScheme, visualizeGridCoords)

        % Method to return the grid coords and the corresponding RTVFobj indices
        % for theConesNumPooled case
        [gridCoords, theRTVFobjIndicesForThisGrid] = subGridSpatialCoordsForConesNumPooled(theConesNumPooled, ...
            theConesNumPooledByTheRFcenterGrid, theOpticsPositionGrid);


        % Method to visualize a computed RTVF object 
        peekIntoSingleRTVFobj(theRTVFTobj, iRTVobjIndex, ...
            theOpticsPositionGrid, theConesNumPooledByTheRFcenterGrid, figNo);


        % Method to visualize the fitted locations
        visualizeFittedLocations(mosaicDirectory, figNo, theMidgetRGCmosaic, ...
            theOpticsPositionGrid, theConesNumPooledByTheRFcenterGrid);


        % Method to visualize the fitted locations
        visualizeFittedLocationsCombo(mosaicDirectory, figNo, theMidgetRGCmosaic, ...
            theRTFVTobjList, theOpticsPositionGrid, ...
            theConesNumPooledByTheRFcenterGrid, varargin);


    end % Class methods
    
end
