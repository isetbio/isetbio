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

        % If not empty, we generate optics at this fixed positions no matter
        % where the RTVF is fitted within the mosaic. If the
        % stfComputeMethod is set to RTVF.directSTFcomputeMethod, which
        % means that the STF is derived from cone mosaic STF responses,
        % this must be set to the eccentricity of the optics under which the
        % cone mosaic STF responses were generated
        fixedZernikeCoefficientsPositionDegs;

        % The nominal spatial sampling grid
        nominalSpatialSamplingGrid;

        % The full sampling grids
        opticalPositionGrid;
        conesNumPooledByTheRFcenterGrid;
        visualSTFSurroundToCenterRcRatioGrid;
        visualSTFSurroundToCenterIntegratedSensitivityRatioGrid;

        % The target RGC indices for each RTVFobj
        targetRGCindex;
        targetRGCindexOfDifferentMajorityConeType;

        % The list with all computed RTVF objects
        RTVFTobjList;

        % Multistarts num
        multiStartsNumRetinalPooling;
        useParallelMultiStart;
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

            % Parse optional input
            p = inputParser;
            p.addParameter('multiStartsNumRetinalPooling', [], @(x)(isempty(x)||isscalar(x)));
            p.addParameter('minSpatialSamplingDegs', 0.25, @isnumeric);
            p.addParameter('useParallelMultiStart', false, @islogical);
            p.addParameter('fixedZernikeCoefficientsPositionDegs', [], @(x)(isempty(x)||((isnumeric(x))&&(numel(x)==2))));
            p.parse(varargin{:});

            obj.multiStartsNumRetinalPooling = p.Results.multiStartsNumRetinalPooling;
            obj.useParallelMultiStart = p.Results.useParallelMultiStart;
            obj.fixedZernikeCoefficientsPositionDegs = p.Results.fixedZernikeCoefficientsPositionDegs;

            if (strcmp(obj.stfComputeMethod, RTVF.directSTFcomputeMethod))
               
                if (isempty(obj.fixedZernikeCoefficientsPositionDegs))
                    error('You must specify the (x,y) eccentricity under which the cone mosaic STF responses were generated')
                end

                % Load the file with the mosaic responses exists
                midgetRGCMosaicInspector.say('Select file with cone mosaic STF responses');

                dropboxDir = midgetRGCMosaicInspector.localDropboxPath();
                [file,path] = uigetfile(fullfile(dropboxDir, '*.mat'), ...
                            'Select a file');
    
                if (file ~= 0)
                    coneMosaicResponsesFileName = fullfile(path,file);
                    fprintf('Will fit RTVFs based on %s\n', coneMosaicResponsesFileName);
                    fprintf(2,'\nAre you sure that the cone mosaic STF responses in this file were generated with optics at position (%2.1f,%2.1f) degs??\n', ...
                            obj.fixedZernikeCoefficientsPositionDegs(1), obj.fixedZernikeCoefficientsPositionDegs(2));
                    fprintf('Hit enter to continue ...')
                    pause;
                    obj.loadComputeMethodResources(coneMosaicResponsesFileName);
                else
                    error('RTVF cannot proceed without cone mosaic STF responses');
                end
            end

            
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

        % Method to determine the majority center cone type for an RGC
        [theCenterConeTypeWeights, theCenterConeTypeNum, theMajorityConeType] = centerConeTypeWeights(obj, theRGCindex);

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