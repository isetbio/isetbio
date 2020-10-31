classdef mRGCmosaic < handle
% Create a midget RGC mosaic connected to a cone mosaic

    % Read-only properties
    properties (GetAccess=public, SetAccess=private)
        % The input cone mosaic
        inputConeMosaic;
    end
    
    % Private properties
    properties (GetAccess=private, SetAccess=private)
        % The metadata of the input cone mosaic (used for connecting cones 
        % to the subregions of the mRGC cell receptive fields.
        inputConeMosaicMetaData;
    end
    
    % Public methods
    methods
        % Constructor
        function obj = mRGCmosaic(eccentricityDegs, sizeDegs, whichEye, varargin)
            % Parse input
            switch (nargin)
                case 0
                    eccentricityDegs = 0;
                    sizeDegs = 1.0;
                    whichEye = 'right';
                case 1
                    sizeDegs = 1.0;
                    whichEye = 'right';
                case 2
                    % do nothing
                    whichEye = 'right';
                otherwise
                    % Parse varargin
            end
            
            % Parse varargin, which could contain an actual cone mosaic
            if (isempty(varargin))
                % An actual cone mosaic was not passed in varargin, so generate one that is appropriate for the eccentricity and size of the mRGC mosaic 
                % Compute cone and mRGC RF positions
                [coneRFpositionsMicrons, coneRFpositionsDegs, rgcRFpositionsMicrons,  rgcRFpositionsDegs, extraDegsForRGCSurround] = ...
                    mRGCmosaic.importConeAndRGCpositions(eccentricityDegs, sizeDegs, whichEye);
               
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
        end
    end
    
    % Static methods
    methods (Static)
        % Static method to import cone and mRGC positions from pre-computed
        % lattices that have the desired size and are centered at the desired eccentricity in the desired eye
        [coneRFpositionsMicrons, coneRFpositionsDegs, ...
         rgcRFpositionsMicrons,  rgcRFpositionsDegs, extraDegsForRGCSurround] = importConeAndRGCpositions(eccentricityDegs, sizeDegs, whichEye);
        
        % Static method to generate a cone mosaic from the imported cone positions
        [theConeMosaic, theConeMosaicMetaData] = generateInputConeMosaic(generationMode, eccentricityDegs, sizeDegs, extraDegsForRGCSurround, coneRFpositionsMicrons);
    end
end

