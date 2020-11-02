classdef mRGCmosaic < handle
% Create a midget RGC mosaic connected to a cone mosaic

    properties (Constant)
        LCONE_ID = 2;
        MCONE_ID = 3;
        SCONE_ID = 4;
    end
    
    
    % Read-only properties
    properties (GetAccess=public, SetAccess=private)
        
        % The eccentricity of the mosaic, in degrees
        eccentricityDegs;
        
        % The size of the mosaic, in degrees
        sizeDegs;
        
        % Eye, left or right
        whichEye;
        
        % The input cone mosaic
        inputConeMosaic;
        
        % Sparse matrix [nCones x mRGC] storing the connection state
        % between the n-th cone to m-th RGC (1==connected, 0==disconencted)
        coneConnectivityMatrix;
        
        % [m x 2] matrix of RGC positions, in microns
        rgcRFpositionsMicrons;
        
        % [m x 2] matrix of RGC positions, in degrees
        rgcRFpositionsDegs;
    end
    
    % Private properties
    properties (GetAccess=private, SetAccess=private)
        % The metadata of the input cone mosaic (used for connecting cones 
        % to the subregions of the mRGC cell receptive fields.
        inputConeMosaicMetaData;
        
        % Size (in degs) of source lattice from which to crop positions for
        % the desired eccentricity
        sourceLatticeSizeDegs = 45;
        
        % [m x 1] matrix of local spacing for each RGC, in microns
        rgcRFspacingsMicrons;
        
        % [m x 1] matrix of local spacing for each RGC, in degs
        rgcRFspacingsDegs
    end
    
    % Public methods
    methods
        % Constructor
        function obj = mRGCmosaic(eccentricityDegs, sizeDegs, whichEye, varargin)
            % Parse input
            switch (nargin)
                case 0
                    eccentricityDegs = 0;
                    sizeDegs = 0.2;
                    whichEye = 'right';
                case 1
                    sizeDegs = 0.5;
                    whichEye = 'right';
                case 2
                    % do nothing
                    whichEye = 'right';
                otherwise
                    % Parse varargin
            end
            
            % Set properties
            obj.eccentricityDegs = eccentricityDegs;
            obj.sizeDegs = sizeDegs;
            obj.whichEye = whichEye;
            
            % 
            % To do: parse varargin, which could contain an actual cone mosaic
            %
            
            if (isempty(varargin))
                % An actual cone mosaic was not passed in varargin, so generate one that is appropriate for the eccentricity and size of the mRGC mosaic 
                % Compute cone and mRGC RF positions
                [coneRFpositionsMicrons, coneRFpositionsDegs, rgcRFpositionsMicrons, rgcRFpositionsDegs, extraDegsForRGCSurround] = ...
                    mRGCmosaic.importConeAndRGCpositions(obj.sourceLatticeSizeDegs, eccentricityDegs, sizeDegs, whichEye);
               
                % Generate a regular hex mosaic to serve as the
                % input cone mosaic with a custom mean cone spacing 
                % (equal to the median spacing within the
                % coneRFpositionsMicrons) and custom quantal efficiency and
                % macular pigment appropriate for the eccentricityDegs
                generationMode = 'equivalent regular hex';
                [obj.inputConeMosaic, obj.inputConeMosaicMetaData] = mRGCmosaic.generateInputConeMosaic(generationMode, ...
                    eccentricityDegs, sizeDegs, extraDegsForRGCSurround, coneRFpositionsMicrons);
                
                plotPositions = true;
                if (plotPositions)
                    coneRFpositionsDegsInRegHexMosaic = RGCmodels.Watson.convert.rhoMMsToDegs(obj.inputConeMosaicMetaData.conePositionsMicrons*1e-3);
                    maxPosDegsX = max(coneRFpositionsDegsInRegHexMosaic(:,1));
                    minPosDegsX = min(coneRFpositionsDegsInRegHexMosaic(:,1));
                    maxPosDegsY = max(coneRFpositionsDegsInRegHexMosaic(:,2));
                    minPosDegsY = min(coneRFpositionsDegsInRegHexMosaic(:,2));
                    figure(1);
                    clf;
                    subplot(1,2,1);
                    plot(coneRFpositionsDegs(:,1), coneRFpositionsDegs(:,2), 'k.');
                    hold on;
                    plot(rgcRFpositionsDegs(:,1), rgcRFpositionsDegs(:,2),'ro');
                    set(gca, 'XLim', [minPosDegsX maxPosDegsX], 'YLim', [minPosDegsY maxPosDegsY]);
                    axis 'equal';
                    title('imported cone positions');
                    
                    subplot(1,2,2);
                   
                    plot(coneRFpositionsDegsInRegHexMosaic(:,1), coneRFpositionsDegsInRegHexMosaic(:,2), 'k.');
                    hold on;
                    plot(rgcRFpositionsDegs(:,1), rgcRFpositionsDegs(:,2),'ro');
                    set(gca, 'XLim', [minPosDegsX maxPosDegsX], 'YLim', [minPosDegsY maxPosDegsY]);
                    axis 'equal';
                    title('regular hex mosaic cone positions');
                end
            end
            
            % Wire cones to RGC center subregions with a cone specificity level
            coneSpecificityLevel = 100;
            
            % Wire cones to RGC centers
            [obj.coneConnectivityMatrix, ...
             obj.rgcRFpositionsDegs, ...
             obj.rgcRFpositionsMicrons, ...
             obj.rgcRFspacingsMicrons] = mRGCmosaic.wireInputConeMosaicToRGCcenters(...
                rgcRFpositionsDegs, rgcRFpositionsMicrons,  ...
                obj.inputConeMosaicMetaData.conePositionsDegs, ...
                obj.inputConeMosaicMetaData.conePositionsMicrons, ...
                obj.inputConeMosaicMetaData.coneSpacingsMicrons, ...
                obj.inputConeMosaicMetaData.coneTypes, ...
                obj.inputConeMosaicMetaData.indicesOfConesNotConnectingToRGCcenters, ...
                coneSpecificityLevel);
            
            % Compute local spacings from positions
            obj.rgcRFspacingsDegs = RGCmodels.Watson.convert.positionsToSpacings(obj.rgcRFpositionsDegs);
        end
        
        % Method to visualize the tesselation of the input cone mosaic by
        % the RF centers of the RGC mosaic
        visualizeConeMosaicTesselation(obj, domain);
    end
    
    % Static methods
    methods (Static)
        % Static method to import cone and mRGC positions from pre-computed
        % lattices that have the desired size and are centered at the desired eccentricity in the desired eye
        [coneRFpositionsMicrons, coneRFpositionsDegs, ...
         rgcRFpositionsMicrons,  rgcRFpositionsDegs, extraDegsForRGCSurround] = ...
            importConeAndRGCpositions(sourceLatticeSizeDegs, eccentricityDegs, sizeDegs, whichEye);
        
        % Static method to generate a cone mosaic from the imported cone positions
        [theConeMosaic, theConeMosaicMetaData] = ...
            generateInputConeMosaic(generationMode, eccentricityDegs, sizeDegs, extraDegsForRGCSurround, coneRFpositionsMicrons); 
        
        % Static method to wire cones to the the RGC RF centers
        [connectivityMatrix, rgcRFpositionsDegs, rgcRFpositionsMicrons, rgcRFspacingsMicrons] = ...
            wireInputConeMosaicToRGCcenters(rgcRFpositionsDegs, rgcRFpositionsMicrons, ...
            conePositionsDegs, conePositionsMicrons, coneSpacingsMicrons, coneTypes, ...
            indicesOfConesNotConnectingToRGCcenters, coneSpecificityLevel);
    end
end

