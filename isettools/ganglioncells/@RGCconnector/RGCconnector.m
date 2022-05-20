classdef RGCconnector < handle
    % A class to handle connection of a lattice of RGCs to a cMosaic
    %
    % Syntax:
    %   theInputConeMosaic = cMosaic();
    %   rc = RGCconnector(theInputConeMosaic);
    %
    %   theConeMosaic = cMosaic('pigment', pp);
    %   theConeMosaic = cMosaic('sizeDegs',[5,5],'eccentricityDegs',[0,0]);
    %
    % Description
    %

     % Public properties
    properties  (GetAccess=public, SetAccess=public)
        
    end % Public properties

    % Read-only properties
    properties (GetAccess=public, SetAccess=private)
        inputConeMosaic;
        RGCRFpositionsMicrons;
        RGCRFspacingsMicrons;
        coneToRGCDensityRatioComputeStruct;
    end

    properties (Constant)
        validRGClatticeModels = {...
            'Watson-midgetRGC' ...
            'Watson-parasolRGC' ...
            };
    end % Constant properties


    % Public methods
    methods
        
        % Constructor
        function obj = RGCconnector(inputConeMosaic, varargin)
            % Validate inputConeMosaic
            assert(isa(inputConeMosaic, 'cMosaic'), 'inputConeMosaic must be a @cMosaic object');
            obj.inputConeMosaic = inputConeMosaic;
            
            % Parse optional input
            p = inputParser;
            p.addParameter('RGCRFpositionsMicrons', [], @(x)((isempty(x)) || (isnumeric(x)&&(size(x,2)==2))));
            p.addParameter('coneToRGCDensityRatio', [], @(x)((isempty(x)) || isnumeric(x)));
            p.parse(varargin{:});
            
            RGCRFposMicrons = p.Results.RGCRFpositionsMicrons;
            coneToRGCDensityRatio = p.Results.coneToRGCDensityRatio;

            if (isempty(RGCRFposMicrons)) && (isempty(coneToRGCDensityRatio))
                % Nothing was specified, so we initialize with a precomputed RGC lattice (Watson's model)
                modelRGC = 'Watson-midgetRGC';
                fprintf('Instantiating using ''%s'' mRGC lattice', modelRGC);
                RGCRFposMicrons = obj.initializeWithPrecomputedLattice(modelRGC);

            elseif (isempty(RGCRFposMicrons)) && (~isempty(coneToRGCDensityRatio))
                % Only a density was specified, so we initialize with a perfect hex RGCRF lattice
                RGCRFposMicrons = obj.initializeWithPerfectHexLattice(coneToRGCDensityRatio);

            elseif (~isempty(RGCRFposMicrons)) && (isempty(coneToRGCDensityRatio))
                % Perhaps issue a warning if there is no overlap with the cMosaic
                fprintf('Instantiating using supplied [%d x 2] lattice of RGC positions\n', ...
                    size(RGCRFposMicrons,1));
            else
                % User specified both a density and a lattice of
                % RGCRFpositions, so raise an error 
                error('@RGCconnector cannot be instantiated both with a set of RGCRFpositions AND a coneToRGCdensityRatio.');
            end

            % Before cropping compute the local cone to RGC density struct
            samplingIntervalMicrons = 3;
            obj.computeConeToRGCDensityRatioComputeStruct(RGCRFposMicrons, samplingIntervalMicrons);

            % Crop positions to lie within the inputConeMosaic
            obj.cropLattice(RGCRFposMicrons);

            % Visualize the input mosaics
            obj.visualizeInputMosaics();

            % Visualize the cone-to-RGC code density map
            obj.visualizeInputConeToRGCDensityMap();

            % Visualize effective lattice and cone to RGC density map
            obj.visualizeEffectiveConeToRGCDensityMap();

        end % Constructor

        % Visualization of input cone mosaics (before any connections are made)
        [hFig, ax, XLims, YLims] = visualizeInputMosaics(obj, varargin)

        % Visualization of cone-to-RGC density map
        visualizeInputConeToRGCDensityMap(obj);

        visualizeEffectiveConeToRGCDensityMap(obj);
        
    end % Public methods


    methods (Access=private)
        % Initialize RGCRF positions by importing them from a large
        % previously-computed mesh based on some RGC model
        RGCRFposMicrons = initializeWithPrecomputedLattice(obj, modelRGC);

        % Initialize RGCRF positions with a perfect hexagonal lattice
        % with a desired coneToRGCDensityRatio
        RGCRFposMicrons = initializeWithPerfectHexLattice(obj, coneToRGCDensityRatio);

        % Crop the imported RGCRFlattice to lie within the inputConeMosaic
        cropLattice(obj, RGCRFposMicrons);

        % Compute the obj.coneToRGCDensityRatioComputeStruct. This is done once, just before cropping.
        computeConeToRGCDensityRatioComputeStruct(obj, RGCRFposMicrons, samplingIntervalMicrons);
    end


    methods (Static)
        [f,v] = facesAndVertices(positions, spacings, diskOutline);
    end
end

