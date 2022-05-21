function visualizeRGCinputs(obj, ax, varargin)
% % Visualize the input cones to each RGC

    % Parse input
    p = inputParser;
    p.addParameter('superimposeConeInputWiring', true, @islogical);
    p.addParameter('displayRGCID', false, @islogical);
    p.addParameter('cMap', brewermap(1024, '*greys'), @isnumeric);
    p.parse(varargin{:});

    superimposeConeInputWiring = p.Results.superimposeConeInputWiring;
    displayRGCID = p.Results.displayRGCID;
    cMap = p.Results.cMap;

    % Find all 0-, 1-, and 2-input RGCs
    zeroInputRGCindices = numel(find(squeeze(sum(obj.coneConnectivityMatrix,1)) == 0))
    oneInputRGCindices = numel(find(squeeze(sum(obj.coneConnectivityMatrix,1)) == 1))
    twoInputRGCindices = numel(find(squeeze(sum(obj.coneConnectivityMatrix,1)) == 2))
    

    rgcsNum = size(obj.coneConnectivityMatrix,2);
    for iRGC = 1:rgcsNum
        % Retrieve connection data
        % Indices of cones connected to this RGC
        connectedConeIndices = find(squeeze(obj.coneConnectivityMatrix(:, iRGC))>0);
        % Weights of these input cones
        inputConeWeights = full(obj.coneConnectivityMatrix(connectedConeIndices, iRGC));

        % Positions and spacings of these input cones
        inputConePositions = obj.inputConeMosaic.coneRFpositionsMicrons(connectedConeIndices,:);
        inputConeSpacings = obj.inputConeMosaic.coneRFspacingsMicrons(connectedConeIndices);

        if (isempty(connectedConeIndices))
            fprintf(2,'Skip rendering RGC %d which has 0 inputs. Drawing an X at its initial position (%2.0f,%2.0f).\n', ...
                iRGC, obj.RGCRFpositionsMicrons(iRGC,1), obj.RGCRFpositionsMicrons(iRGC,2));
            plot(ax, obj.RGCRFpositionsMicrons(iRGC,1), obj.RGCRFpositionsMicrons(iRGC,2), ...
                'h', 'MarkerFaceColor', [0 0 0], 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 16, 'LineWidth', 1);
            continue;
        end

        
        % Generation RF visualization struct
        rfVisualizationStruct = rfFromConeInputs(...
                inputConePositions, inputConeSpacings, inputConeWeights);
        

        % Plot contour of RF
        zLevels(1) = 0.05*min(inputConeWeights);
        zLevels(2) = max(inputConeWeights);
        faceAlpha = 0.2;
        lineStyle = '-';
        lineWidth = 0.5;
        RGCconnector.transparentContourPlot(ax, ...
            rfVisualizationStruct.spatialSupportXY, ...
            rfVisualizationStruct.rfMap, ...
            zLevels, faceAlpha, cMap, lineStyle, lineWidth);


        if (superimposeConeInputWiring)
            if (numel(inputConeWeights)>1)
                centroid = obj.RGCRFcentroidsFromInputs(iRGC,:);
                if (displayRGCID)
                    text(ax, centroid(1), centroid(2), sprintf('%d', iRGC), 'Color', [0 1 0], 'FontSize', 12, 'BackgroundColor', [0 0 0]);
                end
                for iInputCone = 1:numel(inputConeWeights)
                    xx = [centroid(1) inputConePositions(iInputCone,1)];
                    yy = [centroid(2) inputConePositions(iInputCone,2)];
                    plot(ax, xx, yy, 'k-', 'LineWidth', inputConeWeights(iInputCone)/max(inputConeWeights)*3);
                end
            end
        end

    end % iRGC

end


function rfVisualizationStruct = rfFromConeInputs(...
            inputConePositions, inputConeSpacings, inputConeWeights, varargin)
% Generate the 2D RF of an RGC from its input cones
%
% Syntax:
%   [spatialSupportXY, RF, ...
%    theInputConeLineSpreadFunctionsXY] = ...
%         RGCRFconnector.rfFromConeInputs(...
%            inputConePositions, inputConeSpacings, inputConeWeights, varargin)
%
% Description:
%   Generate the 2D RF of an RGC from its input cones
%
% Inputs:
%    inputConePositions   - Vector with indices of the input cones
%    inputConeSpacings    - Vector with the spacing of the input cones
%    inputConeWeights     - Vector with the pooling weights of the input cones
%
%
% Outputs:
%    spatialSupportXY                   - [N x 2] matrix of RF spatial support (x/y)
%    RF                                 - [N x N] matrix of RF
%    theInputConeLineSpreadFunctionsXY  - {2 x inputsNum] matrix of the line spread
%                                         functions of the input cone apertures along 
%                                         the X- and the Y-axes
%
% Optional key/value pairs
%   'xSupport'                  - Custom spatial support (X)
%   'ySupport'                  - Custom spatial support (Y)
%   'spatialSupportSamples'     - Custom spatial support samples
%
% History:
%   5/11/2022       NPC     Wrote it
%

    p = inputParser;
    p.addOptional('xSupport', [], @(x)(isempty(x) || (isnumeric(x))));
    p.addOptional('ySupport', [], @(x)(isempty(x) || (isnumeric(x))));
    p.addOptional('spatialSupportSamples', 60, @(x)(isempty(x) || (isnumeric(x))));
    p.parse(varargin{:});
    xSupport = p.Results.xSupport;
    ySupport = p.Results.ySupport;
    spatialSupportSamples = p.Results.spatialSupportSamples;

    % Compute spatial support
    xSep = max(inputConeSpacings);

    if (isempty(xSupport))
        xx = inputConePositions(:,1);
        xSupport = linspace(min(xx)-xSep,max(xx)+xSep,spatialSupportSamples);
    end

    if (isempty(ySupport))
        yy = inputConePositions(:,2);
        ySupport = linspace(min(yy)-xSep,max(yy)+xSep,spatialSupportSamples);
    end
    [X,Y] = meshgrid(xSupport, ySupport);
    spatialSupportXY(:,1) = xSupport(:);
    spatialSupportXY(:,2) = ySupport(:);
    RF = zeros(size(X));

    % Output
    theInputConeLineSpreadFunctionsXY = cell(2, numel(inputConeWeights));

    for iCone = 1:numel(inputConeWeights)
        % Characteristic radius of the input cone
        rC = 0.204*sqrt(2.0)*inputConeSpacings(iCone);
        % Compute aperture2D x weight
        XX = X-inputConePositions(iCone,1);
        YY = Y-inputConePositions(iCone,2);
        theConeAperture2D = inputConeWeights(iCone) * exp(-(XX/rC).^2) .* exp(-(YY/rC).^2);
        % Accumulate 2D apertures
        RF = RF + theConeAperture2D;
        % 1D Line spread functions (XY) for each cone aperture
        theInputConeLineSpreadFunctionsXY{1,iCone} = sum(theConeAperture2D,1);
        theInputConeLineSpreadFunctionsXY{2,iCone} = sum(theConeAperture2D,2);
    end
    
    % Assemble rfVisualizationStruct
    rfVisualizationStruct.spatialSupportXY = spatialSupportXY;
    rfVisualizationStruct.rfMap = RF;
    rfVisualizationStruct.theInputConeLineSpreadFunctionsXY = theInputConeLineSpreadFunctionsXY;
end
