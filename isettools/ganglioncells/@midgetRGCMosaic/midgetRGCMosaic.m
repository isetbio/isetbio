classdef midgetRGCMosaic < handle
% Create a midgetRGCMosaic mosaic object

    % Public properties
    properties  (GetAccess=public, SetAccess=public)
        % Name of the mosaic
        name;

    end % public properties

    % Read-only properties
    properties (GetAccess=public, SetAccess=private)

        % The source lattice specification
        sourceLatticeSizeDegs;

        % Eccentricity [x,y] of the center of the mRGC mosaic (degs)
        eccentricityDegs;

        % Eccentricity [x,y] of the center of the mRGC mosaic (microns)
        eccentricityMicrons;

        % Size of the mRGC mosaic [width, height]
        sizeDegs;

        % Either 'right eye' or 'left eye'
        whichEye;

        % The input cone mosaic (@cMosaic)
        inputConeMosaic;

        % Extra size of the input cone mosaic (to allow for cone inputs to
        % the RF surrounds)
        extraDegsForInputConeMosaic;

        % Overlap of RFs
        rfOverlapRatio;

        % [0: minimize chromatic variance 1: minimize spatial variance]
        chromaticSpatialVarianceTradeoff = 1.0;

        % Computed params
        % Sparse [conesNum x rgcRFsNum] sparse  connectivity matrix 
        % To find which cones are connected to a n rgcRF:
        %  connectivityVector = full(squeeze(rgcRFcenterConeConnectivityMatrix(:, rgcRF)));
        %  inputConeIndices = find(abs(connectivityVector) > 0.0001);
        rgcRFcenterConeConnectivityMatrix;

        % [n x 2] matrix of rgc rf positions, in microns
        rgcRFpositionsMicrons;
        
        % [n x 1] vector of rgc rf spacings, in microns
        rgcRFspacingsMicrons;

    end % Read-only properties


    % Private properties
    properties (GetAccess=private, SetAccess=private)
    
        
    end


    % Public methods
    methods
        
        % Constructor
        function obj = midgetRGCMosaic(varargin)
            % Parse input
            p = inputParser;
            p.addParameter('name', 'mRGC mosaic', @ischar);
            p.addParameter('inputConeMosaic', [], @(x)(isempty(x) || isa(x, 'cMosaic')));
           
            p.addParameter('sourceLatticeSizeDegs', 58, @isscalar);
            p.addParameter('whichEye', 'right eye', @(x)(ischar(x) && (ismember(x, {'left eye', 'right eye'}))));
            
            % Eccentricity/Position: if positionDegs is sent in, it overwrites eccentricityDegs
            p.addParameter('eccentricityDegs', [], @(x)(isnumeric(x) && (numel(x) == 2)));
            p.addParameter('positionDegs', [], @(x)(isnumeric(x) && (numel(x) == 2)));
            p.addParameter('sizeDegs', [0.4 0.4], @(x)(isnumeric(x) && (numel(x) == 2)));

            % Custom Degs <-> MM conversion functions
            p.addParameter('customDegsToMMsConversionFunction', @(x)RGCmodels.Watson.convert.rhoDegsToMMs(x), @(x) (isempty(x) || isa(x,'function_handle')));
            p.addParameter('customMMsToDegsConversionFunction', @(x)RGCmodels.Watson.convert.rhoMMsToDegs(x), @(x) (isempty(x) || isa(x,'function_handle')));
            
            % RF overlap
            p.addParameter('rfOverlapRatio', 0, @(x)(isscalar(x)&&((x>=0)&&(x<=1))));
            p.addParameter('chromaticSpatialVarianceTradeoff', 1.0, @(x)(isscalar(x)&&((x>=0)&&(x<=1))));

            p.parse(varargin{:});

            obj.name = p.Results.name;
            obj.sourceLatticeSizeDegs = p.Results.sourceLatticeSizeDegs;
            obj.rfOverlapRatio = p.Results.rfOverlapRatio;
            obj.chromaticSpatialVarianceTradeoff = p.Results.chromaticSpatialVarianceTradeoff;


            % Deal with the input cone mosaic. If an input cMosaic is not specified we generate one here
            if (~isempty(p.Results.inputConeMosaic))
                fprintf(2,'An input cone mosaic was passed: ignoring set key/value pairs for properties inherited by the cMosaic.\n')
                pause
            else
                % Generate the input cone mosaic. This also sets various
                % other properties that are matched to the input cone
                % mosaic, such as:
                % - eccentricityDegs
                % - sizeDegs
                % - whichEye
                obj.generateInputConeMosaic(p.Results);
                
                % Set the eccentricity in microns
                obj.eccentricityMicrons = obj.inputConeMosaic.eccentricityMicrons;

                % Wire the RF centers. This sets the following properties:
                % - rgcRFcenterConeConnectivityMatrix
                % - rgcRFpositionsMicrons
                % - rgcRFspacingsMicrons
                obj.generateRFpositionsAndWireTheirCenters();

            end
        end % Constructor
    end % Public methods

    % Private methods
    methods (Access=private)
        generateInputConeMosaic(obj, pResults);
        generateRFpositionsAndWireTheirCenters(obj);
    end % Private methods

    % Static methods
    methods (Static)
        extraDegs = extraConeMosaicDegsForMidgetRGCSurrounds(eccentricityDegs, sizeDegs);
    end % Static metthods

end
